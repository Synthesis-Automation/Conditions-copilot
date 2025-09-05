from __future__ import annotations
from typing import Dict, Any, List

def _find_base_name(reagent_order: List[Dict[str, Any]] | None) -> str:
    if not reagent_order:
        return "?"
    # Prefer explicit role match
    for r in reagent_order:
        if isinstance(r, dict) and r.get("role") == "base":
            name = r.get("name")
            if isinstance(name, str) and name:
                return name
    # Fallback: first item that looks like a base by common patterns
    for r in reagent_order:
        name = r.get("name") if isinstance(r, dict) else None
        if isinstance(name, str) and any(t in name for t in ("CO3", "OtBu", "PO4")):
            return name
    return "?"

def predict_yield_stub(reaction_smiles: str, conditions: Dict[str, Any]) -> Dict[str, float]:
    """Deterministic, defensive stub returning an uncalibrated score.
    - Looks at `condition_core` and presence of carbonate/alkoxide base.
    - Never raises on missing or short reagent lists.
    """
    core = conditions.get("condition_core") or "?"
    fc = conditions.get("full_conditions") or {}
    order = fc.get("reagent_order") if isinstance(fc, dict) else None
    base = _find_base_name(order if isinstance(order, list) else [])

    bias = 0.5
    if isinstance(core, str):
        if "Pd" in core:
            bias = 0.7
        elif "Ni" in core:
            bias = 0.6
    if isinstance(base, str) and "CO3" in base:
        bias += 0.05

    pred = round(float(bias), 2)
    succ = min(0.95, max(0.3, float(bias)))
    return {"predicted_yield": pred, "success_prob": succ}
