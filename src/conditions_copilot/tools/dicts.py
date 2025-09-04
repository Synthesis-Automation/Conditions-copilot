from __future__ import annotations
from typing import Dict, Any
import os, json

DEFAULTS = {
    "catalysts": ["Pd2(dba)3","Pd(OAc)2","NiCl2(dme)","CuI"],
    "ligands": ["SPhos","XPhos","BrettPhos","dtbbpy","1,10-phenanthroline"],
    "bases": ["Cs2CO3","K3PO4","NaOtBu","K2CO3"],
    "solvents": ["1,4-dioxane","toluene","DMA","DMSO"],
    "additives": ["DMAP","TBAB"]
}

def load_dicts(dir_path: str) -> Dict[str, Any]:
    out = {}
    for key in DEFAULTS.keys():
        p = os.path.join(dir_path, f"{key}.json")
        if os.path.exists(p):
            out[key] = json.load(open(p))
        else:
            out[key] = DEFAULTS[key]
    return out
