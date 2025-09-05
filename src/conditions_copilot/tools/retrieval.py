from __future__ import annotations
from typing import List, Dict, Any, Tuple, Set
import os

import pandas as pd

# Minimal neighbor stub: demo neighbors for Ar–Br + aniline if dataset missing.
DEMO_NEIGHBORS = [
    {"ConditionCore": "Pd/SPhos", "tails": "Cs2CO3/1,4-dioxane/90-100°C", "yield": 0.78,
     "substrate_tags": ["Ar–Br", "aniline"], "notes": "Robust for Ar–Br + aniline"},
    {"ConditionCore": "Pd/XPhos", "tails": "K3PO4/toluene/100-110°C", "yield": 0.75,
     "substrate_tags": ["Ar–Br", "aniline"], "notes": "Bulky ligand; add 10-20% DMA if needed"},
    {"ConditionCore": "Ni/dtbbpy", "tails": "NaOtBu/1,4-dioxane/70-80°C", "yield": 0.68,
     "substrate_tags": ["Ar–Br", "aniline"], "notes": "Stronger base"},
]


def _parse_bits(s: str) -> Set[int]:
    if not isinstance(s, str) or not s:
        return set()
    out: Set[int] = set()
    for tok in s.split("|"):
        tok = tok.strip()
        if not tok:
            continue
        try:
            out.add(int(tok))
        except Exception:
            continue
    return out


def _tanimoto(a: Set[int], b: Set[int]) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return float(inter) / float(union) if union else 0.0


def retrieve_neighbors(dataset_csv: str, k: int = 25, filters: Dict[str, Any] | None = None
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Neighbor retrieval with optional ECFP4-diff Tanimoto ranking.

    If the dataset CSV has a column `ecfp4_diff_bits` ("|"-separated bit indices),
    and `filters` provides `query_bits`, we rank by Tanimoto similarity to the
    query and break ties with yield. Otherwise we fall back to yield sorting.
    """
    if not os.path.exists(dataset_csv):
        return DEMO_NEIGHBORS[:min(k, len(DEMO_NEIGHBORS))], {
            "n_total": 0, "n_close": len(DEMO_NEIGHBORS), "chemotype_overlap": "high", "median_yield_close": 0.74
        }

    df = pd.read_csv(dataset_csv)

    sub = df.copy()
    # Optional tag filter pass (compatible with existing CSVs)
    if filters and "class_tags" in filters and "substrate_tags" in df.columns:
        tags = set(filters["class_tags"]) if isinstance(filters["class_tags"], (list, set, tuple)) else {filters["class_tags"]}

        def match(row):
            st = set(str(row.get("substrate_tags", "")).split("|"))
            return bool(st & tags)

        sub = sub[sub.apply(match, axis=1)].copy()

    # Tanimoto path
    sims = None
    if "ecfp4_diff_bits" in df.columns and filters and "query_bits" in filters:
        try:
            qbits: Set[int] = {int(i) for i in filters["query_bits"] if str(i).strip()}
            sims = sub["ecfp4_diff_bits"].fillna("").map(lambda s: _tanimoto(qbits, _parse_bits(s)))
        except Exception:
            sims = None

    if sims is not None:
        sub = sub.assign(_sim=sims).sort_values(["_sim", "yield"], ascending=[False, False]).head(k)
    else:
        sub = sub.sort_values("yield", ascending=False).head(k)

    neighbors: List[Dict[str, Any]] = []
    for _, row in sub.iterrows():
        neighbors.append({
            "ConditionCore": row.get("ConditionCore", "unknown"),
            "tails": row.get("tails", ""),
            "yield": float(row.get("yield", 0.0)),
            "substrate_tags": str(row.get("substrate_tags", "Ar–X|amine")).split("|"),
            "notes": row.get("notes", ""),
        })

    cov = {
        "n_total": int(len(df)),
        "n_close": int(len(sub)),
        "chemotype_overlap": "high" if "ecfp4_diff_bits" in df.columns else "medium",
        "median_yield_close": float(sub["yield"].median()) if len(sub) else None,
    }
    return neighbors, cov

