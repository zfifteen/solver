# Milestone 15 Release Candidate

This document records the acceptance surface for the current release-candidate state. For the fuller capability inventory, limitations, and milestone snapshot, see [PROJECT_STATUS.md](./PROJECT_STATUS.md).

## Supported Release Surface

- Apple Silicon macOS only
- deterministic CPU-first solver path as the reference validation contract
- restartable lid-driven cavity simulation path via the versioned checkpoint format
- targeted Metal support for 3D periodic Taylor-Green only

## Acceptance Sources

| Area | Source of Truth |
| --- | --- |
| Project status and supported capabilities | [PROJECT_STATUS.md](./PROJECT_STATUS.md) |
| Deterministic operating procedure | [DETERMINISTIC_RUNBOOK.md](./DETERMINISTIC_RUNBOOK.md) |
| Validation and benchmark evidence | [validation/latest/summary.md](../validation/latest/summary.md), [validation/latest/benchmark_summary.md](../validation/latest/benchmark_summary.md) |
| Performance evidence | [profiling/latest/summary.md](../profiling/latest/summary.md), [profiling/latest/throughput_summary.md](../profiling/latest/throughput_summary.md), [profiling/latest/hotspot_summary.md](../profiling/latest/hotspot_summary.md) |
| Reproducible suite entrypoints | [tools/run_release_candidate_suite.sh](../tools/run_release_candidate_suite.sh), `ctest --test-dir build/deterministic --output-on-failure`, [validation/run_validation_suite.py](../validation/run_validation_suite.py) |

## Current Acceptance Snapshot

- Validation overall status is `PASS`.
- Benchmark thresholds and mass-conservation gates are passing on the current reference evidence set.
- The current Apple M1 Max profiling snapshot recommends the `benchmark` build with the default unclamped scheduler policy.
- The current compute-thread recommendation remains `1` while the solver path is still effectively single-threaded.
- The supported restartable full-simulation path is the lid-driven cavity driver.
- The targeted Metal backend is validated only for the supported 3D periodic Taylor-Green slice.

## Final Acceptance Checklist

| Item | Status |
| --- | --- |
| Root documentation reflects the release-candidate surface | Pass |
| Deterministic configure, build, test, and validation flow is documented | Pass |
| Fresh validation artifacts under `validation/latest/` pass | Pass |
| Fresh profiling artifacts under `profiling/latest/` pass | Pass |
| Restartable cavity path is documented and covered by automated tests | Pass |

## Acceptance Rule

Milestone 15 is satisfied when the repository state, the deterministic runbook, and the generated validation and profiling artifacts agree on the supported release surface.
