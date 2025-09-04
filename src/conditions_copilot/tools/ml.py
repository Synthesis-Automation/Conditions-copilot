from __future__ import annotations
from typing import Dict, Any

def predict_yield_stub(reaction_smiles: str, conditions: Dict[str, Any]) -> Dict[str, float]:
    """Tiny deterministic stub: returns an uncalibrated score just for wiring tests."""
    core = conditions.get("condition_core","?")
    base = (conditions.get("full_conditions",{}).get("reagent_order",[]) or [{}])[2].get("name","?") if conditions.get("full_conditions") else "?"
    bias = 0.7 if "Pd" in core else 0.6 if "Ni" in core else 0.5
    if "CO3" in base: bias += 0.05
    return {"predicted_yield": round(bias, 2), "success_prob": min(0.95, max(0.3, bias))}
