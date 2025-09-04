# Conditions Copilot (Repo Skeleton)

A Copilot-style, multi-iteration **reaction condition recommender**. The local app is a thin
bootstrapper + tool runner. The remote LLM does discovery/planning and calls your tools.

## Features
- **Universal discovery payload** (reaction-agnostic).
- **Typed JSON schemas** (pydantic v2).
- **Tool stubs**: `featurize_basic`, `retrieve_neighbors`, `get_dicts`, `predict_yield`, `validate_proposals`.
- **Validator** (L0–L3): schema/units → dictionaries/constraints → compatibility/solubility → dataset support.
- **Controller loop** with budgets (≤3 tool calls per cycle; ≤5 cycles).
- **CLI** for discovery payload and basic iteration (interactive by default).

## Quickstart
```bash
python -m venv .venv && source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# 1) Print a universal discovery payload for a reaction
python -m conditions_copilot.cli discover --rxn "Brc1ccccc1C(=O)C.Nc1ccccc1>>CC(=O)c1ccccc1Nc1ccccc1"

# 2) Run an interactive iteration (paste LLM JSON replies when prompted)
python -m conditions_copilot.cli run --rxn "Brc1ccccc1C(=O)C.Nc1ccccc1>>CC(=O)c1ccccc1Nc1ccccc1"
```

> To integrate an external LLM CLI (e.g., Codex), set `LLM_CMD` environment variable to a shell command that reads stdin and outputs JSON. The controller will pipe the payload into it.

## Layout
```
conditions-copilot/
  ├─ src/conditions_copilot/
  │   ├─ capabilities.py         # Tool manifest for discovery
  │   ├─ schemas.py              # Pydantic models for envelopes & proposals
  │   ├─ controller.py           # Iteration controller (discovery → propose → validate loop)
  │   ├─ cli.py                  # Simple CLI
  │   ├─ llm_prompting/
  │   │   ├─ system.txt          # System rules (planner/critic; JSON-only)
  │   │   └─ client.py           # LLM adapter (env-driven; interactive fallback)
  │   └─ tools/
  │       ├─ featurize_basic.py  # RDKit-lite featurizer (works even if RDKit missing)
  │       ├─ retrieval.py        # kNN stub (CSV-based); demo neighbors if dataset absent
  │       ├─ dicts.py            # Allow-list loaders (with defaults)
  │       ├─ ml.py               # Tiny yield predictor stub
  │       └─ validator.py        # L0–L3 validator + ValidationReport
  ├─ data/dicts/*.json           # Default dictionaries
  ├─ examples/                   # Example payload/expected replies
  ├─ requirements.txt
  └─ README.md
```

## Notes
- RDKit is optional here; if not present, `featurize_basic` returns minimal signals and the loop still works.
- Replace demo neighbors with your dataset (`CSV/Parquet`) in `tools/retrieval.py`.
- Tight schemas + validator keep the LLM honest (JSON-only, allow-lists, safety rules).
