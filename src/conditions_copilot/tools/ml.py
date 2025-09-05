from __future__ import annotations
from typing import Dict, Any, List, Optional
import os, pickle

# Optional deps
try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
except Exception:
    Chem = None
    DataStructs = None
    AllChem = None

try:
    from sklearn.base import BaseEstimator  # type: ignore
except Exception:
    BaseEstimator = object  # type: ignore

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


def _ecfp4_diff_bits_from_rxn(rxn_smiles: str, n_bits: int = 2048) -> Optional[List[int]]:
    if Chem is None or AllChem is None:
        return None
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
        if bv_r is None or bv_p is None or DataStructs is None:
            return None
        bv_d = (bv_p ^ bv_r)
        return list(bv_d.GetOnBits())
    except Exception:
        return None


def _vectorize(reaction_smiles: str, conditions: Dict[str, Any], n_bits: int = 2048) -> Dict[str, Any]:
    """Lightweight featurizer for ML inference.
    Returns a dict with:
      - fp_bits: list[int] (ECFP4 diff bit indices)
      - knobs: simple numeric knobs from conditions if present (temp, time)
    """
    bits = _ecfp4_diff_bits_from_rxn(reaction_smiles, n_bits=n_bits) or []
    knobs = {}
    try:
        fc = conditions.get("full_conditions", {}) if isinstance(conditions, dict) else {}
        if isinstance(fc, dict):
            if "temperature_c" in fc:
                knobs["temperature_c"] = float(fc.get("temperature_c"))
            if "time_h" in fc:
                knobs["time_h"] = float(fc.get("time_h"))
    except Exception:
        pass
    return {"fp_bits": bits, "knobs": knobs}


def _to_dense(vec: Dict[str, Any], n_bits: int = 2048) -> List[float]:
    arr = [0.0] * (n_bits + 2)
    for i in vec.get("fp_bits", []):
        try:
            ii = int(i)
        except Exception:
            continue
        if 0 <= ii < n_bits:
            arr[ii] = 1.0
    # Append simple knobs in last slots
    t = vec.get("knobs", {}).get("temperature_c")
    h = vec.get("knobs", {}).get("time_h")
    arr[-2] = float(t) if isinstance(t, (int, float)) else 0.0
    arr[-1] = float(h) if isinstance(h, (int, float)) else 0.0
    return arr


def predict_yield(reaction_smiles: str, conditions: Dict[str, Any]) -> Dict[str, float]:
    """Enhanced predictor (optional). Falls back cleanly if model/deps missing.

    Looks for a scikit-learn regressor pickle at path from env `YIELD_MODEL_PATH`
    or `models/yield_model.pkl` relative to the repo. Expects `.predict` and
    optionally `.predict_proba`. Feature vector: 2048-bit ECFP4 diff + temp/time.
    """
    model_path = os.environ.get("YIELD_MODEL_PATH") or os.path.join(os.path.dirname(__file__), "..", "..", "models", "yield_model.pkl")
    try:
        if not os.path.exists(model_path):
            raise FileNotFoundError
        with open(model_path, "rb") as f:
            model = pickle.load(f)
        vec = _vectorize(reaction_smiles, conditions)
        x = [_to_dense(vec)]
        yhat = None
        succ = None
        if hasattr(model, "predict"):
            try:
                yhat = float(model.predict(x)[0])
            except Exception:
                yhat = None
        if hasattr(model, "predict_proba"):
            try:
                proba = model.predict_proba(x)[0]
                if isinstance(proba, (list, tuple)):
                    succ = float(proba[-1])
                else:
                    succ = float(proba)
            except Exception:
                succ = None
        if yhat is None and succ is not None:
            yhat = max(0.0, min(1.0, succ * 0.9))
        if yhat is None:
            raise RuntimeError("model_predict_failed")
        return {"predicted_yield": round(float(yhat), 2), "success_prob": round(float(succ if succ is not None else yhat), 2)}
    except Exception:
        # Fallback to rule-based stub
        return predict_yield_stub(reaction_smiles, conditions)
