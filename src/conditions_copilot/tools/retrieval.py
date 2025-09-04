from __future__ import annotations
from typing import List, Dict, Any, Tuple
import os, json, pandas as pd

# Minimal neighbor stub: demo neighbors for Ar–Br + aniline if dataset missing.
DEMO_NEIGHBORS = [
    {"ConditionCore":"Pd/SPhos","tails":"Cs2CO3/1,4-dioxane/90–100°C","yield":0.78,
     "substrate_tags":["Ar–Br","aniline"],"notes":"Robust for Ar–Br + aniline"},
    {"ConditionCore":"Pd/XPhos","tails":"K3PO4/toluene/100–110°C","yield":0.75,
     "substrate_tags":["Ar–Br","aniline"],"notes":"Bulky ligand; add 10–20% DMA if needed"},
    {"ConditionCore":"Ni/dtbbpy","tails":"NaOtBu/1,4-dioxane/70–90°C","yield":0.68,
     "substrate_tags":["Ar–Br","aniline"],"notes":"Stronger base"},
]

def retrieve_neighbors(dataset_csv: str, k: int = 25, filters: Dict[str, Any] | None = None
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    if not os.path.exists(dataset_csv):
        return DEMO_NEIGHBORS[:min(k, len(DEMO_NEIGHBORS))], {
            "n_total": 0, "n_close": len(DEMO_NEIGHBORS), "chemotype_overlap": "high", "median_yield_close": 0.74
        }
    df = pd.read_csv(dataset_csv)
    # Very naive retrieval: filter by class_tags if present; then take top-k by yield
    if filters and "class_tags" in filters:
        tags = set(filters["class_tags"])
        def match(row):
            st = set(str(row.get("substrate_tags","" )).split("|"))
            return bool(st & tags)
        sub = df[df.apply(match, axis=1)].copy()
    else:
        sub = df.copy()
    sub = sub.sort_values("yield", ascending=False).head(k)
    neighbors = []
    for _, row in sub.iterrows():
        neighbors.append({
            "ConditionCore": row.get("ConditionCore","unknown"),
            "tails": row.get("tails",""),
            "yield": float(row.get("yield", 0.0)),
            "substrate_tags": str(row.get("substrate_tags","Ar–X|amine")).split("|"),
            "notes": row.get("notes",""),
        })
    cov = {"n_total": int(len(df)), "n_close": int(len(sub)),
           "chemotype_overlap": "medium", "median_yield_close": float(sub["yield"].median()) if len(sub) else None}
    return neighbors, cov
