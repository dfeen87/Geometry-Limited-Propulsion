# Verification Plan (Step-Up for External Review)

This repository already defines falsifiable claims. This plan makes verification easier to execute and audit.

## Scope

- Validate numerical behavior of the gate and Δv integration code.
- Validate reproducibility of notebook-derived claims.
- Prepare minimum artifacts for independent review.

## Level 1 — Code Reliability (immediate)

1. Create isolated environment.
2. Install package in editable mode.
3. Run unit tests.
4. Record Python/package versions and commit hash.

**Pass criteria**
- All tests pass.
- No import-path workarounds required.
- Results are repeatable on two consecutive runs.

## Level 2 — Numerical Consistency (short-term)

1. Run deterministic sweeps for `coherence_gate` and `hysteresis_sweep`.
2. Recompute selected notebook figures from script/notebook.
3. Confirm monotonic/threshold behavior matches README claims.

**Pass criteria**
- Curves and summary metrics remain stable under fixed seed/config.
- No sign or unit inconsistencies.

## Level 3 — Experimental Readiness (medium-term)

1. Pre-register perturbation protocol and metrics.
2. Define sensor/calibration requirements for Δφ proxies.
3. Publish raw data schema + analysis script interface.

**Pass criteria**
- Independent reviewer can run analysis from raw data to headline plots.
- Falsification conditions are mechanically checkable.

## Recommended Artifacts to Add Next

- `verification/baseline_results.json` (expected summary outputs)
- `verification/reproduce_notebook.py` (headless reproduction)
- `verification/data_schema.md` (required columns, units, sampling rates)

## Minimum Evidence Bundle for Aerospace Review

- Commit hash + environment manifest
- Test log
- Figure/table reproduction log
- Explicit list of assumptions and unvalidated mappings
