# Conditions Copilot — Phased Implementation Plan

This roadmap breaks improvements into concrete phases with goals, key tasks, deliverables, and acceptance criteria. It reflects work already completed and what’s next, optimized for iterative shipping and measurable gains.

## Scope & Assumptions
- Python app with CLI and Qt GUI, local tools + optional LLM planner.
- RDKit optional; graceful degradation required.
- Datasets are CSV-first; may include precomputed fingerprints.

---

## Phase 0 — Baseline Stabilization (Complete / Ongoing)
- Goals: Improve UX, fix brittle JSON parsing, keep app usable without RDKit.
- Done:
  - GUI: colorized transcript, dark theme, higher contrast.
  - Lenient JSON repair for LLM replies (code fences, smart quotes, trailing commas, Python-literal fallback).
  - Enhanced featurizer behind optional RDKit; schema now allows extra feature keys.
- Deliverables:
  - Resilient GUI and parser; richer features shipped to LLM.
- Acceptance:
  - GUI renders colored roles + readable content on dark background.
  - Paste/external LLM replies parse successfully in typical failure modes.

---

## Phase 1 — Feature Engineering Upgrade (Shipped, Extend)
- Goals: Stronger substrate/reaction features with RDKit, graceful fallback without RDKit.
- Shipped:
  - ECFP4 unions and reaction-diff bits (2048) and descriptor aggregates + deltas.
  - Rule-oriented features inspired by Rule-conditions: amine, carboxylic acid derivatives, aryl halides.
- Next tasks:
  - Add reaction-center SMARTS for common classes (amide, ester, Suzuki, Buchwald–Hartwig).
  - Protonation flags (pH 8–10) and “likely salt” booleans.
- Acceptance:
  - New fields present in features; no crashes when RDKit unavailable.

---

## Phase 2 — Retrieval Upgrade (Shipped, Improve)
- Goals: Similarity-aware neighbors + robust coverage metrics.
- Shipped:
  - Tanimoto ranking using `ecfp4_diff_bits` when present; yield-sort fallback.
  - Precompute utility: `python -m conditions_copilot.tools.precompute_fps input.csv`.
- Next tasks:
  - Auto-wire live `ECFP4_diff_bits` from featurizer into retrieval as `filters.query_bits` when the LLM asks for neighbors.
  - Diversify top-k (cluster or MaxMin) and expose similarity quantiles in coverage.
- Acceptance:
  - Retrieval ranks by Tanimoto when bits provided; otherwise matches prior behavior.

---

## Phase 3 — ML Predictor (Shipped Fallback, Add Model)
- Goals: Replace rule stub with a trained, calibrated scorer; safe fallback retained.
- Shipped:
  - `predict_yield` loads `models/yield_model.pkl` (or `YIELD_MODEL_PATH`) and falls back to `predict_yield_stub`.
  - Feature vector: 2048-bit ECFP4 diff + temp/time knobs.
- Next tasks:
  - Add `scripts/train_yield_model.py` (RF/XGBoost), calibration (Platt/Isotonic), and export pipeline.
  - Save model and version metadata; doc install instructions.
- Acceptance:
  - Predictor returns `predicted_yield` and `success_prob` with sensible calibration on holdout.

---

## Phase 4 — LLM Loop: Tool Schemas + Critic Retry
- Goals: Make the loop more reliable and self-correcting.
- Tasks:
  - Add explicit tool JSON Schemas to prompts; validate replies; compact error diff to LLM on failure; one auto-retry within budget.
  - Optional few-shot with 2–3 nearest neighbor exemplars (when coverage high).
- Acceptance:
  - Reduced invalid-proposal rate; measurable increase in first-pass valid responses.

---

## Phase 5 — Panel Design (DoE)
- Goals: Generate robust 2×3 or 3×3 screening panels within dictionaries/constraints.
- Tasks:
  - New `tools/panel.py` for factorial or D-optimal selection over ligand/base/solvent/temperature.
  - Include diversity penalty to avoid near-duplicates.
- Acceptance:
  - Panel fields present when coverage low or uncertainty high.

---

## Phase 6 — Data Quality & Synonyms
- Goals: Better normalization and negative data support.
- Tasks:
  - Expand `data/dicts` with aliases and canonical names; add normalizer in `tools/dicts.py`.
  - Ingest failed/low-yield cases; label for calibration.
- Acceptance:
  - Proposals strictly use canonical dict items; model calibration improves.

---

## Phase 7 — Observability & Reproducibility
- Goals: Inspectability and consistent reruns.
- Tasks:
  - Save run artifacts to `runs/<timestamp>/` (inputs, tool outputs, LLM replies, validation).
  - Add simple metrics (validator pass rate, predictor MAE/ROC if labels available).
- Acceptance:
  - Runs are reproducible; artifacts aid debugging and iteration.

---

## Phase 8 — GUI Enhancements
- Goals: Decision support in the UI.
- Tasks:
  - Inline predictions for proposed conditions + uncertainty badges; validator reasons as tooltips.
  - Theme toggle (dark/light) and quick filters.
- Acceptance:
  - Users can compare proposals with predicted yield/risks at a glance.

---

## Phase 9 — Performance & Caching
- Goals: Faster, cheaper, repeatable.
- Tasks:
  - Tool result caching keyed by `(tool, args_hash)`; optional disk cache with TTL.
  - Batch operations where possible (e.g., predict multiple proposals in one call).
- Acceptance:
  - Re-running typical flows avoids recompute; latency decreases noticeably.

---

## Phase 10 — Safety & Validation Tightening
- Goals: Prevent unsafe or out-of-scope suggestions.
- Tasks:
  - Expand validator checks (solvent/base incompatibilities, temperature windows, metal blacklists, safety flags).
  - Stop-ship blocker thresholds; structured feedback to LLM.
- Acceptance:
  - Validator catches unsafe/infeasible proposals consistently; LLM self-corrects on retry.

---

## Milestones & Suggested Timeline
- Week 1–2: Phase 2 auto-wire `query_bits`, Phase 4 schemas + critic retry, Phase 9 caching.
- Week 3–4: Phase 3 training script + calibrated model, Phase 5 DoE panels.
- Week 5–6: Phase 6 synonyms + negatives ingestion, Phase 7 observability, Phase 8 GUI inline scores.

---

## Risks & Mitigations
- RDKit availability: keep optional; gate features; provide clear install docs.
- Dataset sparsity/bias: blend NN ex with model priors; expose uncertainty; suggest panels.
- LLM drift: schema-validated tool calls; critic retries; minimal deterministic prompting.

---

## Next Actions (Checklist)
- [ ] Auto-inject `ECFP4_diff_bits` into retrieval filters (`controller.execute_actions`).
- [ ] Add tool result cache (in-memory + optional disk).
- [ ] Add tool schemas + critic retry in `llm_prompting/client.py`.
- [ ] Implement `tools/panel.py` and surface panel in proposals.
- [ ] Add `scripts/train_yield_model.py` and model export docs.
- [ ] Save run artifacts to `runs/<timestamp>/` and minimal metrics.

References (Where to edit)
- Featurization: `src/conditions_copilot/tools/featurize_basic.py`
- Retrieval: `src/conditions_copilot/tools/retrieval.py`
- ML: `src/conditions_copilot/tools/ml.py`, `models/`
- Controller: `src/conditions_copilot/controller.py`
- GUI: `src/conditions_copilot/gui_qt.py`
- LLM client: `src/conditions_copilot/llm_prompting/client.py`

