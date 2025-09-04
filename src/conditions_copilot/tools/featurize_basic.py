from __future__ import annotations
from typing import Dict, Any
import math

# RDKit is optional; we degrade gracefully
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except Exception:
    Chem = None
    Descriptors = None

def featurize_basic(smiles: str) -> Dict[str, Any]:
    """Return coarse, reaction-agnostic features.
    Works even if RDKit is missing (minimal signals).
    """
    lhs, _, rhs = smiles.partition(">>")
    reactants = [s for s in lhs.split(".") if s]
    feats: Dict[str, Any] = {
        "HALIDE_CLASS": None,
        "NUCLEOPHILE_CLASS": None,
        "BORON_CLASS": None,
        "OXIDATION_DELTA": False,
        "REDUCTION_DELTA": False,
        "C_C_FORMATION": ("C.C" in lhs) and ("C" in rhs),  # silly placeholder
        "C_N_FORMATION": ("N" in lhs and "N" in rhs),
        "ORTHO_COUNT": None,
        "EWG_FLAGS": [],
        "RING_STATS": {"n_aromatic_rings_reactants": None, "n_aromatic_rings_products": None},
        "TPSA": None,
        "logP": None,
    }
    if Chem is None:
        return feats  # minimal signals only

    def mols(ss): return [Chem.MolFromSmiles(s) for s in ss if s]
    r_mols = mols(reactants)
    p_mols = mols([s for s in rhs.split(".") if s])

    # simple HALIDE_CLASS
    hal = None
    for m in r_mols:
        for a in m.GetAtoms():
            if a.GetSymbol() in ("Cl","Br","I") and any(n.GetIsAromatic() for n in a.GetNeighbors()):
                hal = {"Cl":"Ar–Cl","Br":"Ar–Br","I":"Ar–I"}[a.GetSymbol()]; break
        if hal: break
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
        if nuc: break
    feats["NUCLEOPHILE_CLASS"] = nuc

    # TPSA/logP
    if r_mols:
        feats["TPSA"] = round(sum(Descriptors.TPSA(m) for m in r_mols), 2)
        feats["logP"] = round(sum(Descriptors.MolLogP(m) for m in r_mols), 2)

    # ring stats (rough)
    def arom_count(ms):
        c = 0
        for m in ms:
            for ring in m.GetRingInfo().AtomRings():
                if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    c += 1
        return c
    feats["RING_STATS"]["n_aromatic_rings_reactants"] = arom_count(r_mols) if r_mols else None
    feats["RING_STATS"]["n_aromatic_rings_products"] = arom_count(p_mols) if p_mols else None

    return feats
