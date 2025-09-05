from __future__ import annotations
from typing import Dict, Any, List
import math

# RDKit is optional; we degrade gracefully
try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors
except Exception:
    Chem = None
    DataStructs = None
    Descriptors = None
    AllChem = None
    rdMolDescriptors = None


def _mols(smiles_list: List[str]):
    if Chem is None:
        return []
    ms = []
    for s in smiles_list:
        try:
            m = Chem.MolFromSmiles(s)
            if m is not None:
                ms.append(m)
        except Exception:
            continue
    return ms


def _bv_union(mols, radius=2, n_bits=2048):
    if AllChem is None:
        return None
    acc = None
    for m in mols:
        try:
            bv = AllChem.GetMorganFingerprintAsBitVect(m, radius=radius, nBits=n_bits)
        except Exception:
            continue
        acc = bv if acc is None else (acc | bv)
    return acc


def _onbits(bv) -> List[int] | None:
    if bv is None:
        return None
    try:
        return list(bv.GetOnBits())
    except Exception:
        return None


def featurize_basic(smiles: str) -> Dict[str, Any]:
    """Return reaction-agnostic features with optional RDKit extras.
    - Always returns minimal signals even without RDKit.
    - With RDKit: adds ECFP4 unions/diff and richer descriptors incl. deltas.
    """
    lhs, _, rhs = smiles.partition(">>")
    reactants = [s for s in lhs.split(".") if s]
    products = [s for s in rhs.split(".") if s]

    feats: Dict[str, Any] = {
        # coarse placeholders (no-RDKit compatible)
        "HALIDE_CLASS": None,
        "NUCLEOPHILE_CLASS": None,
        "BORON_CLASS": None,
        "OXIDATION_DELTA": False,
        "REDUCTION_DELTA": False,
        "C_C_FORMATION": ("C.C" in lhs) and ("C" in rhs),
        "C_N_FORMATION": ("N" in lhs and "N" in rhs),
        "ORTHO_COUNT": None,
        "EWG_FLAGS": [],
        "RING_STATS": {"n_aromatic_rings_reactants": None, "n_aromatic_rings_products": None},
        "TPSA": None,
        "logP": None,
    }

    if Chem is None:
        return feats  # minimal signals only

    r_mols = _mols(reactants)
    p_mols = _mols(products)

    # simple HALIDE_CLASS
    hal = None
    for m in r_mols:
        for a in m.GetAtoms():
            if a.GetSymbol() in ("Cl", "Br", "I") and any(n.GetIsAromatic() for n in a.GetNeighbors()):
                hal = {"Cl": "Ar–Cl", "Br": "Ar–Br", "I": "Ar–I"}[a.GetSymbol()]
                break
        if hal:
            break
    feats["HALIDE_CLASS"] = hal

    # crude amine class
    nuc = None
    for m in r_mols:
        for a in m.GetAtoms():
            if a.GetSymbol() == "N":
                if any(nb.GetIsAromatic() for nb in a.GetNeighbors()):
                    nuc = "aniline (primary aromatic amine)"
                else:
                    nuc = "amine"
                break
        if nuc:
            break
    feats["NUCLEOPHILE_CLASS"] = nuc

    # TPSA/logP
    if r_mols:
        try:
            feats["TPSA"] = round(sum(Descriptors.TPSA(m) for m in r_mols), 2)
            feats["logP"] = round(sum(Descriptors.MolLogP(m) for m in r_mols), 2)
        except Exception:
            pass

    # ring stats (rough)
    def arom_count(ms):
        c = 0
        for m in ms:
            try:
                for ring in m.GetRingInfo().AtomRings():
                    if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                        c += 1
            except Exception:
                continue
        return c

    feats["RING_STATS"]["n_aromatic_rings_reactants"] = arom_count(r_mols) if r_mols else None
    feats["RING_STATS"]["n_aromatic_rings_products"] = arom_count(p_mols) if p_mols else None

    # Rich descriptors (aggregates and deltas)
    def agg_desc(ms):
        d = {
            "MolWt": 0.0,
            "HBA": 0,
            "HBD": 0,
            "RotB": 0,
            "RingCount": 0,
            "FracCSP3": 0.0,
            "AromFrac": 0.0,
        }
        if not ms:
            return d
        n_atoms = 0
        arom_atoms = 0
        for m in ms:
            try:
                d["MolWt"] += float(Descriptors.MolWt(m))
                d["HBA"] += int(rdMolDescriptors.CalcNumHBA(m))
                d["HBD"] += int(rdMolDescriptors.CalcNumHBD(m))
                d["RotB"] += int(Descriptors.NumRotatableBonds(m))
                d["RingCount"] += int(Descriptors.RingCount(m))
                d["FracCSP3"] += float(Descriptors.FractionCSP3(m))
                n_atoms += m.GetNumAtoms()
                arom_atoms += sum(1 for a in m.GetAtoms() if a.GetIsAromatic())
            except Exception:
                continue
        d["FracCSP3"] = round(d["FracCSP3"], 3)
        d["AromFrac"] = round((arom_atoms / n_atoms) if n_atoms else 0.0, 3)
        d["MolWt"] = round(d["MolWt"], 2)
        return d

    r_desc = agg_desc(r_mols)
    p_desc = agg_desc(p_mols)
    delta = {k + "_DELTA": round(p_desc[k] - r_desc[k], 3) for k in r_desc}
    feats["DESCRIPTORS_REACTANTS"] = r_desc
    feats["DESCRIPTORS_PRODUCTS"] = p_desc
    feats["DESCRIPTORS_DELTA"] = delta

    # ECFP4 unions and diff
    bv_r = _bv_union(r_mols, radius=2, n_bits=2048)
    bv_p = _bv_union(p_mols, radius=2, n_bits=2048)
    bv_diff = None
    if bv_r is not None and bv_p is not None and DataStructs is not None:
        try:
            bv_diff = (bv_p ^ bv_r)
        except Exception:
            bv_diff = None
    feats["ECFP4_reactants_bits"] = _onbits(bv_r)
    feats["ECFP4_products_bits"] = _onbits(bv_p)
    feats["ECFP4_diff_bits"] = _onbits(bv_diff)

    return feats
