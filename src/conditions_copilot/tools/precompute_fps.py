from __future__ import annotations
import argparse, os
from typing import List

import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except Exception:
    Chem = None
    AllChem = None


def ecfp4_diff_bits_from_rxn(rxn_smiles: str, n_bits: int = 2048) -> List[int]:
    if Chem is None or AllChem is None:
        return []
    lhs, _, rhs = rxn_smiles.partition(">>")
    r_ss = [s for s in lhs.split(".") if s]
    p_ss = [s for s in rhs.split(".") if s]
    try:
        def bv_union(ss):
            acc = None
            for s in ss:
                m = Chem.MolFromSmiles(s)
                if m is None:
                    continue
                bv = AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=n_bits)
                acc = bv if acc is None else (acc | bv)
            return acc
        bv_r = bv_union(r_ss)
        bv_p = bv_union(p_ss)
        if bv_r is None or bv_p is None:
            return []
        bv_d = (bv_p ^ bv_r)
        return list(bv_d.GetOnBits())
    except Exception:
        return []


def main():
    ap = argparse.ArgumentParser(description="Precompute ECFP4 diff bit indices for reactions in a CSV.")
    ap.add_argument("input_csv", help="Path to input CSV containing a 'rxn_smiles' column")
    ap.add_argument("output_csv", nargs="?", help="Path to output CSV (default: overwrite input)")
    args = ap.parse_args()

    if Chem is None or AllChem is None:
        raise SystemExit("RDKit not available. Install rdkit to use this utility.")

    inp = args.input_csv
    out = args.output_csv or args.input_csv
    df = pd.read_csv(inp)
    if "rxn_smiles" not in df.columns:
        raise SystemExit("Expected a column 'rxn_smiles' in the dataset.")
    bits_col = []
    for s in df["rxn_smiles"].fillna(""):
        bits = ecfp4_diff_bits_from_rxn(str(s))
        bits_col.append("|".join(str(i) for i in bits))
    df["ecfp4_diff_bits"] = bits_col
    df.to_csv(out, index=False)
    print(f"Wrote {out} with ecfp4_diff_bits")


if __name__ == "__main__":
    main()

