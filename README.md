# Conditions Copilot (Repo Skeleton)

A Copilot-style, multi-iteration reaction condition recommender. The local app is a thin
bootstrapper + tool runner. The remote LLM does discovery/planning and calls your tools.

## Features
- Universal discovery payload (reaction-agnostic).
- Typed JSON schemas (Pydantic v2).
- Tool stubs: `featurize_basic`, `retrieve_neighbors`, `get_dicts`, `predict_yield`, `validate_proposals`.
- Validator (L0–L3): schema/units → dictionaries/constraints → compatibility/solubility → dataset support.
- Controller loop with budgets (≤3 actions per cycle; ≤4 cycles).
- CLI for discovery payload and basic iteration (interactive by default).

## Installation

Option A — Editable install (recommended):
```bash
python -m venv .venv && source .venv/bin/activate   # Windows (PowerShell): .venv\Scripts\Activate.ps1
pip install -e .
```

Option B — Without install (use PYTHONPATH):
```bash
python -m venv .venv && source .venv/bin/activate   # Windows (PowerShell): .venv\Scripts\Activate.ps1
pip install -r requirements.txt

# macOS/Linux
export PYTHONPATH=src
# Windows (PowerShell)
$env:PYTHONPATH = 'src'
```

## Quickstart
```bash
# 1) Print a universal discovery payload for a reaction
python -m conditions_copilot.cli discover --rxn "Brc1ccccc1C(=O)C.Nc1ccccc1>>CC(=O)c1ccccc1Nc1ccccc1"

# 2) Run an interactive iteration (paste LLM JSON replies when prompted)
python -m conditions_copilot.cli run --rxn "Brc1ccccc1C(=O)C.Nc1ccccc1>>CC(=O)c1ccccc1Nc1ccccc1"
```

To integrate an external LLM CLI (e.g., Codex), set `LLM_CMD` to a shell command
that reads stdin and outputs JSON. The controller pipes the payload into it.

## Layout
```
conditions-copilot/
  ├─ src/conditions_copilot/
  │  ├─ capabilities.py         # Tool manifest for discovery
  │  ├─ schemas.py              # Pydantic models for envelopes & proposals
  │  ├─ controller.py           # Iteration controller (discovery → propose → validate)
  │  ├─ cli.py                  # Simple CLI
  │  ├─ llm_prompting/
  │  │  ├─ system.txt          # System rules (planner/critic; JSON-only)
  │  │  └─ client.py           # LLM adapter (env-driven; interactive fallback)
  │  └─ tools/
  │      ├─ featurize_basic.py  # RDKit-lite featurizer (works even if RDKit missing)
  │      ├─ retrieval.py        # kNN stub (CSV-based); demo neighbors if dataset absent
  │      ├─ dicts.py            # Allow-list loaders (with defaults)
  │      ├─ ml.py               # Tiny yield predictor stub
  │      └─ validator.py        # L0–L3 validator + ValidationReport
  ├─ data/dicts/*.json           # Default dictionaries
  ├─ examples/                   # Example payload/expected replies and datasets
  │  ├─ datasets/demo_neighbors.csv
  ├─ requirements.txt
  └─ README.md
```

## Notes
- RDKit is optional; if not present, `featurize_basic` returns minimal signals and the loop still works.
- Replace demo neighbors with your dataset (CSV/Parquet) in `tools/retrieval.py`.
- Tight schemas + validator keep the LLM honest (JSON-only, allow-lists, safety rules).

## GUI (Qt6)
- Launch: `conditions-copilot-gui`
- Layout: top transcript (system/tools/LLM responses), bottom input for your reaction SMILES.
- LLM: set `LLM_CMD` or use the GUI “LLM Cmd…” button. If empty, the app prompts you to paste the JSON replies.

## OpenAI-Compatible Integration
- The app talks to OpenAI-compatible chat completion APIs (OpenAI or Aliyun DashScope).
- Set either:
  - `OPENAI_API_KEY` (+ optional `OPENAI_MODEL`, `OPENAI_BASE_URL`), or
  - `DASHSCOPE_API_KEY` (+ optional `DASHSCOPE_MODEL`, `DASHSCOPE_BASE_URL`).
- Defaults if unset:
  - `MODEL`: `deepseek-v3.1`
  - `BASE_URL`: DashScope compatible endpoint when `DASHSCOPE_API_KEY` is present.
- Safety prompt: both CLI and GUI ask for confirmation before sending the payload to the API.
- JSON-only: the client requests JSON responses and enforces JSON parsing.

PowerShell example (DashScope / DeepSeek):

```
$env:DASHSCOPE_API_KEY = '<your key>'
$env:DASHSCOPE_MODEL   = 'deepseek-v3.1'
# Optional (uses this by default when DASHSCOPE_API_KEY is set):
$env:DASHSCOPE_BASE_URL = 'https://dashscope.aliyuncs.com/compatible-mode/v1'
```

## Datasets
- A small example dataset is provided at `examples/datasets/demo_neighbors.csv`.
- CLI default `--dataset` points there; you can override with your own CSV.
