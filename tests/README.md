# Test Suite

A layered safety net so you can validate the pipeline in **seconds/minutes**
instead of rerunning the full multi-day job. None of the 9 numbered scripts are
modified — the suite is purely additive.

## Layers

| Layer | Marker | Speed | What it catches |
|-------|--------|-------|-----------------|
| **Smoke** | `smoke` | seconds | Syntax / import / config-wiring / unresolved-path breakage in every script |
| **Unit** | `unit` | seconds | Logic regressions in the bug-prone pure functions (R/S counting, BLAST tiers, √p colsample, bootstrap CI, threshold, chunking, registry, path resolution) |
| **Integration** | `integration`, `slow` | minutes | The real 01→07 chain on a tiny synthetic dataset (wiring, end-to-end correctness, data-leakage guard) |

## How to run

```bash
# Fast feedback — default. Runs unit + smoke, NEVER touches your config.
pytest

# Just one layer
pytest -m unit
pytest -m smoke

# Full end-to-end on a synthetic 16-genome dataset (organism = "testorg").
# Needs xgboost + optuna + the bundled KMC binary (bin/bin/kmc); otherwise skips.
pytest -m integration -s
```

`pytest` defaults to `-m "not integration and not slow"` (see `pytest.ini`), so
the heavy / config-mutating test is **opt-in only**.

## What the integration test does

1. Generates ~16 deterministic synthetic genomes (8 carry a planted
   "resistance motif") + a binary phenotype CSV — see `conftest.py`.
2. Stages them at the `testorg` organism's resolved paths.
3. Temporarily swaps in a small test config (`k=13`, `chunk_size=4`,
   `n_trials=2`, …) — **restored in a `finally` block**, byte-for-byte.
4. Runs `02 → 07` as subprocesses and asserts each step's key outputs.
5. Includes a **data-leakage guard**: asserts step 06 does NOT overwrite the
   training-derived threshold.
6. Cleans up every `testorg` artifact afterwards.

Steps `08`/`09` (BLAST / Nextflow / NCBI network) are out of scope for the
automated run.

## Environments

In an environment without `xgboost`/`optuna`/`KMC` (e.g. CI lint box), the
xgboost-dependent unit tests and the integration test **skip gracefully**; smoke
+ the remaining unit tests still run. In your real training environment
everything runs.
