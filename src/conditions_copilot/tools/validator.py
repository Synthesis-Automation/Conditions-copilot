from __future__ import annotations
from typing import Dict, Any, List
from . import dicts as dicts_mod

# L0–L3 validator: schema checks are handled by pydantic models upstream.
# Here we do dictionary membership, simple compatibility, and support scoring.

def validate_proposals(proposals_resp: Dict[str, Any], dict_dir: str) -> Dict[str, Any]:
    dicts = dicts_mod.load_dicts(dict_dir)
    flags: List[str] = []
    proposals = proposals_resp.get("proposals", [])
    if not proposals:
        return {"rules_passed": False, "confidence": 0.0, "flags": ["no_proposals"], "support": {}}

    # Dictionary membership checks
    def in_list(name, lst): return (name in lst) if name and isinstance(lst, list) else False

    for p in proposals:
        fc = p.get("full_conditions", {})
        for r in fc.get("reagent_order", []):
            role, name = r.get("role"), r.get("name")
            if role == "catalyst" and not in_list(name, dicts.get("catalysts", [])):
                flags.append(f"item_not_in_dictionary:catalyst:{name}")
            if role == "ligand" and not in_list(name, dicts.get("ligands", [])):
                flags.append(f"item_not_in_dictionary:ligand:{name}")
            if role == "base" and not in_list(name, dicts.get("bases", [])):
                flags.append(f"item_not_in_dictionary:base:{name}")
            if role == "solvent" and not in_list(name, dicts.get("solvents", [])):
                flags.append(f"item_not_in_dictionary:solvent:{name}")

        # Simple base/solvent compatibility warnings
        base_names = [r.get("name") for r in fc.get("reagent_order", []) if r.get("role") == "base"]
        solv_names = [r.get("name") for r in fc.get("reagent_order", []) if r.get("role") == "solvent"]
        if base_names and solv_names:
            b, s = base_names[0], solv_names[0]
            if ("CO3" in b) and (s == "toluene"):
                flags.append("solvent/base_incompat: carbonate in neat toluene (slurry note required)")

    # Support/Confidence aggregation (very rough)
    support = proposals[0].get("support", {})
    neighbor_count = int(support.get("neighbor_count", 0))
    med_y = float(support.get("neighbor_median_yield", 0.0))
    conf = min(1.0, 0.4 + 0.02*neighbor_count + 0.4*med_y)
    rules_passed = not any(f.startswith("item_not_in_dictionary") for f in flags)
    return {
        "rules_passed": rules_passed,
        "confidence": round(conf, 2),
        "flags": flags,
        "support": {"neighbor_count": neighbor_count, "median_yield_close": med_y}
    }
