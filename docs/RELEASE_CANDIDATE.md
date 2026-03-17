# Milestone 15 Release Candidate

This document records the current release-candidate scope, evidence sources, and acceptance state for the solver repository.

## Supported Scope

- Apple Silicon macOS only
- Deterministic CPU-first solver path as the reference validation contract
- Restartable lid-driven cavity simulation path via the versioned checkpoint format
- Narrow Metal support for 3D periodic Taylor-Green only

## Non-Goals

- Broader GPU generalization beyond the current Taylor-Green Metal slice
- Broader 3D full-case coverage for cavity or channel drivers
- Cross-platform portability, packaging, or release tagging
- Expansion beyond the current deterministic benchmark catalog

## Milestone 15 Deliverables

| Deliverable | Source of Truth |
| --- | --- |
| Complete documentation | [`README.md`](../README.md), [`docs/TECH-SPEC.md`](./TECH-SPEC.md), [`docs/EXECUTION_ROADMAP_V1.md`](./EXECUTION_ROADMAP_V1.md), [`docs/DETERMINISTIC_RUNBOOK.md`](./DETERMINISTIC_RUNBOOK.md) |
| Verified benchmark results | [`validation/latest/summary.md`](../validation/latest/summary.md), [`validation/latest/benchmark_summary.md`](../validation/latest/benchmark_summary.md) |
| Performance report | [`profiling/latest/summary.md`](../profiling/latest/summary.md), [`profiling/latest/throughput_summary.md`](../profiling/latest/throughput_summary.md), [`profiling/latest/hotspot_summary.md`](../profiling/latest/hotspot_summary.md) |
| Reproducible test suite | [`tools/run_release_candidate_suite.sh`](../tools/run_release_candidate_suite.sh), `ctest --preset deterministic`, [`validation/run_validation_suite.py`](../validation/run_validation_suite.py) |
| Restartable simulations | [`README.md`](../README.md), [`docs/DETERMINISTIC_RUNBOOK.md`](./DETERMINISTIC_RUNBOOK.md), runtime checkpoint/restart coverage in [`tests/test_runtime.cpp`](../tests/test_runtime.cpp) |
| Documented deterministic build and run procedure | [`docs/DETERMINISTIC_RUNBOOK.md`](./DETERMINISTIC_RUNBOOK.md) |

## Canonical Evidence

### Validation / Benchmark Results

- Overall validation summary: [`validation/latest/summary.md`](../validation/latest/summary.md)
- Benchmark thresholds: [`validation/latest/benchmark_summary.md`](../validation/latest/benchmark_summary.md)
- Supporting tables and plots: `validation/latest/`

Current validation snapshot:

- Overall status: `PASS`
- Operator spatial order gate: `pass`
- Taylor-Green temporal self-convergence gate: `pass`
- Benchmark thresholds gate: `pass`
- Mass-conservation gate: `pass`
- Benchmark highlights:
  - `couette_128`: `profile_relative_l2_error=2.459100e-03`
  - `poiseuille_128`: `profile_relative_l2_error=2.239400e-05`
  - `lid_driven_cavity_re100_128`: `max_reference_relative_error=2.776500e-03`
  - `taylor_green_128`: `normalized_energy_error=7.017980e-04`
  - `taylor_green_3d_64`: `normalized_energy_error=2.488760e-03`

### Performance Report

- Performance summary: [`profiling/latest/summary.md`](../profiling/latest/summary.md)
- Throughput summary: [`profiling/latest/throughput_summary.md`](../profiling/latest/throughput_summary.md)
- Hotspot summary: [`profiling/latest/hotspot_summary.md`](../profiling/latest/hotspot_summary.md)
- Supporting tables and plots: `profiling/latest/`

Current performance snapshot:

- Recommended execution mode: `benchmark` build with the default unclamped scheduler policy
- Compute-thread recommendation: `1`
- Advection microbenchmark throughput: `1.474020e+07` cell updates/s
- Pressure microbenchmark throughput: `1.929090e+06` unknown updates/s
- End-to-end Taylor-Green throughput: `1.101220e+06` cells/s
- Targeted Metal 3D validation result: `1.389x` speedup on `taylor_green_3d_64` with velocity relative L2 drift `1.473779e-03`
- Targeted Metal 256^3 smoke result: `2.208x` speedup on `taylor_green_3d_256_smoke`

## Restart Contract

The supported restartable full-simulation path for this release candidate is the lid-driven cavity driver.

Checkpoint files are:

- versioned
- little-endian
- checksum-verified
- validated against build and configuration hashes on load

The manual restart smoke procedure is documented in [`docs/DETERMINISTIC_RUNBOOK.md`](./DETERMINISTIC_RUNBOOK.md). The automated runtime tests remain the main proof that split cavity runs reproduce uninterrupted deterministic runs and that corrupted or incompatible checkpoints are rejected.

## Final Acceptance Checklist

| Item | Status |
| --- | --- |
| Documentation updated for Milestone 15 release-candidate state | Pass |
| Deterministic configure/build/test flow is reproducibly documented | Pass |
| Deterministic wrapper suite exists for the supported validation surface | Pass |
| Fresh deterministic validation artifacts in `validation/latest/` pass | Pass |
| Fresh profiling artifacts in `profiling/latest/` pass | Pass |
| Restartable cavity path is documented and covered by automated tests | Pass |

## Current Release-Candidate State

Milestone 15 is satisfied when the current repository state, the deterministic runbook, and the generated evidence artifacts all agree. The linked validation and profiling summaries are the canonical numeric record for this release-candidate pass.
