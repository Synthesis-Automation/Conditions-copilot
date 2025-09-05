# Conditions Copilot — Session Notes

This document summarizes key changes and how to run the app, including OpenAI setup and troubleshooting. Useful if you restart your environment.

## What’s New
- Qt6 GUI: `src/conditions_copilot/gui_qt.py`, console entry `conditions-copilot-gui`.
- OpenAI integration: automatic LLM calls when `OPENAI_API_KEY` is set, with a confirm-before-send prompt.
- Example dictionaries: `examples/dicts_example/`.
- Example dataset: `examples/datasets/demo_neighbors.csv`.
- Safer stubs and schema fix: hardened ML stub; `Retrieval.negatives` uses `default_factory=list`.

## Install
- Editable install: `pip install -e .`
- Or: `pip install -r requirements.txt`
- Optional: RDKit extra `pip install -e .[rdkit]`

## Run
- CLI discovery: `python -m conditions_copilot.cli discover --rxn "<SMILES>"`
- CLI loop: `python -m conditions_copilot.cli run --rxn "<SMILES>" --dataset examples/datasets/demo_neighbors.csv --dictdir examples/dicts_example`
- GUI: `conditions-copilot-gui`

## LLM selection order
1. `LLM_CMD` (shell command)
2. `OPENAI_API_KEY` (OpenAI Chat Completions, asks to confirm before sending)
3. Paste mode (manual JSON paste)

To use OpenAI, ensure `LLM_CMD` is unset/empty and `OPENAI_API_KEY` is set.

## OpenAI Setup
- Set: PowerShell `$env:OPENAI_API_KEY='sk-…'`; macOS/Linux `export OPENAI_API_KEY='sk-…'`
- Optional: `OPENAI_MODEL` (default `gpt-4o-mini`), `OPENAI_BASE_URL` (custom endpoint)
- Restart the shell/app after adding env vars via Windows GUI so the process inherits them.

## Troubleshooting
- Check vars: PowerShell `echo $env:OPENAI_API_KEY`; `echo $env:LLM_CMD`
- Confirm dependency: `pip show openai`
- Quick env check: `python -c "import os;print(bool(os.getenv('OPENAI_API_KEY')))`
- Try diag script: `python examples/openai_diag.py`
- Clear `LLM_CMD` if set (it has priority): PowerShell `Remove-Item Env:LLM_CMD` or GUI “LLM Cmd…” → clear
- Network/proxy: ensure outbound HTTPS allowed to OpenAI or your custom base URL
- Model name: set `OPENAI_MODEL` to a valid model (e.g., `gpt-4o-mini`)

## Paths
- Dictionaries: `examples/dicts_example/` (JSON allow-lists)
- Dataset: `examples/datasets/demo_neighbors.csv`

## Notes
- All LLM responses must be JSON; the client enforces JSON parsing.
- GUI shows system payloads, tool results, LLM replies, and validation report in the transcript.
