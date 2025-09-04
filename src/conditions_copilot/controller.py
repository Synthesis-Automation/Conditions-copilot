from __future__ import annotations
import json, os
from typing import Dict, Any
from .capabilities import CAPABILITIES
from .schemas import DiscoveryPayload, ProposalsRequest, ProposalsResponse, Constraints, Reaction, Features, Retrieval, Dictionaries
from .tools import featurize_basic, retrieval, dicts as dicts_mod, validator
from .llm_prompting.client import call_llm

DEFAULT_CONSTRAINTS = Constraints()

def build_discovery_payload(rxn_smiles: str, include_dicts: bool = False) -> Dict[str, Any]:
    payload = {
        "reaction": {"smiles": rxn_smiles, "type_hint": None},
        "context": {"plugins_detected": [], "plugin_confidence": 0.0},
        "capabilities": CAPABILITIES,
        "constraints": DEFAULT_CONSTRAINTS.model_dump(),
        "state": {"cycle": 1, "budget": {"max_steps": 3, "used": 0}},
        "ask": ("1) Classify the reaction (name + confidence). 2) List minimal features and neighbor filters you need. "
                 "3) Emit an ACTION PLAN as JSON (<=3 tool calls). If ambiguous, include a typed yes/no question. Return JSON only.")
    }
    if include_dicts:
        payload["dictionaries"] = dicts_mod.load_dicts(os.path.join(os.path.dirname(__file__), "..", "..", "data", "dicts"))
    return DiscoveryPayload(**payload).model_dump()

def execute_actions(actions: list[dict], rxn_smiles: str, dataset_csv: str, dict_dir: str) -> dict:
    results = {}
    for act in actions[:3]:  # budget
        tool = act.get("tool")
        args = act.get("args", {})
        cid = act.get("call_id") or tool
        if tool == "featurize_basic":
            feats = featurize_basic.featurize_basic(rxn_smiles)
            results[cid] = {"ok": True, "features": feats}
        elif tool == "retrieve_neighbors":
            k = int(args.get("k", 20))
            filters = args.get("filters", {})
            neighbors, coverage = retrieval.retrieve_neighbors(dataset_csv=dataset_csv, k=k, filters=filters)
            results[cid] = {"ok": True, "neighbors": neighbors, "coverage": coverage}
        elif tool == "get_dicts":
            scope = args.get("reaction_class")  # not used to filter in this stub
            dicts = dicts_mod.load_dicts(dict_dir)
            results[cid] = {"ok": True, "dictionaries": dicts}
        elif tool == "predict_yield":
            from .tools.ml import predict_yield_stub
            results[cid] = {"ok": True, **predict_yield_stub(rxn_smiles, args.get("conditions", {}))}
        elif tool == "validate_proposals":
            results[cid] = {"ok": True, "note": "Use local validator entrypoint from CLI."}
        else:
            results[cid] = {"ok": False, "error": f"unknown_tool:{tool}"}
    return results

def run_once(rxn_smiles: str, dataset_csv: str, dict_dir: str, system_path: str) -> dict:
    # Round 1: discovery
    disc = build_discovery_payload(rxn_smiles)
    reply1 = call_llm(disc, system_path)
    actions = reply1.get("actions", [])
    tool_results = execute_actions(actions, rxn_smiles, dataset_csv, dict_dir)

    # Round 2: synthesis/proposals
    # Construct a proposals request from tool_results
    # Find features, neighbors/coverage, dictionaries
    feats = {}
    retr = {"neighbors": [], "coverage": {"n_total":0,"n_close":0,"chemotype_overlap":"low","median_yield_close":None}, "negatives": []}
    dicts = dicts_mod.load_dicts(dict_dir)

    for k, res in tool_results.items():
        if res.get("ok") and "features" in res:
            feats = res["features"]
        if res.get("ok") and "neighbors" in res:
            retr["neighbors"] = res["neighbors"]
            retr["coverage"] = res["coverage"]
        if res.get("ok") and "dictionaries" in res:
            dicts = res["dictionaries"]

    req = {
        "reaction": {"smiles": rxn_smiles, "type_hint": reply1.get("classification",{}).get("name"), "class_confidence": reply1.get("classification",{}).get("confidence")},
        "features": feats,
        "retrieval": retr,
        "dictionaries": dicts,
        "constraints": DEFAULT_CONSTRAINTS.model_dump(),
        "ask": ("Propose top 3 ConditionCore and full conditions; include 2×3 screening panel if weak support. "
                 "Use only dictionary items. Return JSON only.")
    }
    req = ProposalsRequest(**req).model_dump(by_alias=True)
    reply2 = call_llm(req, system_path)

    # Validate
    report = validator.validate_proposals(reply2, dict_dir=dict_dir)
    return {"round1": reply1, "tool_results": tool_results, "request2": req, "round2": reply2, "validation": report}
