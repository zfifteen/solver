# Project Status

This document is the canonical home for the repository's current phase, implemented capabilities, active limitations, and milestone snapshot.

## Current Phase

The repository is currently in **Milestone 15: Release Candidate** for a deliberately narrow supported surface:

- Apple Silicon macOS only
- deterministic CPU-first validation as the reference contract
- restartable lid-driven cavity path
- targeted Metal support for the 3D periodic Taylor-Green path

This is not a general-purpose CFD platform yet. It is a focused solver program that already includes a complete numerical path for its current benchmark surface, plus the validation and profiling machinery needed to support that surface with evidence.

## Supported Scope

- platform: macOS on Apple Silicon with Apple Clang and CMake
- physics: incompressible, constant-density flow
- mesh model: structured Cartesian grids on a staggered MAC layout
- reference execution mode: deterministic CPU-first workflow
- benchmark cases: lid-driven cavity at `Re = 100`, analytic Couette and Poiseuille channel cases, and 2D / 3D Taylor-Green vortex decay
- GPU scope: targeted Metal backend for 3D periodic Taylor-Green only

## Implemented Capabilities

### Platform and Build

- Apple-Silicon-only CMake configuration with labeled `deterministic` and `benchmark` presets
- Apple Clang based build using C++20 plus Objective-C++ for the Metal host path
- smoke executable, minimal CTest wiring, and a release-candidate wrapper script for the supported validation surface

### Numerical Core

- structured `Grid` representation plus MAC-aware pressure, scalar, and velocity storage
- flat contiguous `double` buffers with ghost-cell-aware indexing helpers
- second-order discrete `gradient`, `divergence`, and `laplacian` operators
- transport kernels for advection and diffusion
- bounded TVD advection with `van Leer`, plus first-order upwind and centered alternatives for controlled studies
- pressure-projection machinery with pressure boundary-condition mapping and null-space handling
- semi-implicit momentum path with deterministic ADI line solves
- matrix-free pressure Poisson operator and MGPCG solve with fixed V-cycle geometric multigrid preconditioning

### Simulation and Verification Paths

- deterministic lid-driven cavity timestepper and benchmark driver
- generalized physical-boundary abstraction shared across solver paths
- channel-flow verification driver covering Couette and Poiseuille analytic cases
- Taylor-Green verification driver covering both 2D and 3D periodic cases
- manufactured-solution operator verification tool for reusable spatial convergence studies

### Restart, Output, and Diagnostics

- versioned little-endian binary checkpoint format with checksum verification
- build and configuration hash validation during checkpoint load
- bitwise deterministic cavity continuation and restart coverage
- legacy VTK export for cell-centered pressure and reconstructed velocity
- residual-history tracking, solver diagnostics, and benchmark metadata reporting

### Validation and Profiling Infrastructure

- automated validation harness writing convergence tables, benchmark summaries, conservation summaries, and plots under `validation/latest/`
- automated profiling harness writing kernel, throughput, policy-study, thread-scaling, and hotspot reports under `profiling/latest/`
- Instruments-backed hotspot summaries and comparison reporting against the stored Milestone 10 baseline

### Targeted Metal Backend

- dedicated `metal/` module with Objective-C++ host orchestration and compiled Metal shaders
- explicit `backend=cpu|metal` Taylor-Green configuration and `--backend` CLI override
- comparison reporting in the profiling harness for supported 3D Taylor-Green cases

## Active Limitations

- The cavity and channel full-case drivers remain intentionally 2D-only.
- The cavity and channel drivers remain CPU-only.
- The Metal backend is intentionally narrow: 3D, fully periodic, Taylor-Green-only.
- The GPU path uses float working storage on the current machine because Metal double-precision kernels are not available there; the run finishes with CPU-side cleanup to restore the strict divergence contract.
- The solver path is still effectively single-threaded, so the current thread-scaling study is baseline reporting rather than a mature scaling story.
- Broader 3D full-case coverage, higher-Reynolds-number validation, release packaging, version tagging, and cross-platform distribution are not implemented yet.

## Evidence and Source of Truth

Current operational evidence lives in the generated summaries:

- validation summary: [validation/latest/summary.md](../validation/latest/summary.md)
- benchmark thresholds: [validation/latest/benchmark_summary.md](../validation/latest/benchmark_summary.md)
- profiling summary: [profiling/latest/summary.md](../profiling/latest/summary.md)
- throughput summary: [profiling/latest/throughput_summary.md](../profiling/latest/throughput_summary.md)
- hotspot summary: [profiling/latest/hotspot_summary.md](../profiling/latest/hotspot_summary.md)

Current headline snapshot:

- validation overall status: `PASS`
- operator spatial order gate: `pass`
- Taylor-Green temporal self-convergence gate: `pass`
- benchmark thresholds gate: `pass`
- mass-conservation gate: `pass`
- measured recommendation on the current Apple M1 Max reference machine: `benchmark` build with the default unclamped scheduler policy
- current compute-thread recommendation: `1`
- targeted Metal comparison result: `1.389x` speedup on `taylor_green_3d_64`
- targeted Metal smoke result: `2.208x` speedup on `taylor_green_3d_256_smoke`

For the deterministic operating procedure, see [DETERMINISTIC_RUNBOOK.md](./DETERMINISTIC_RUNBOOK.md). For the release-candidate acceptance snapshot, see [RELEASE_CANDIDATE.md](./RELEASE_CANDIDATE.md).

## Milestone Snapshot

| Milestone Range | Theme | State |
| --- | --- | --- |
| 0-5 | environment, grid/field infrastructure, operators, transport, projection, linear solve | complete |
| 6-8 | first full simulation path, generalized boundary conditions, restart/output | complete |
| 9-11 | broader verification, profiling, and first CPU optimization pass | complete |
| 12-13 | 3D Taylor-Green extension and targeted Metal acceleration | complete |
| 14-15 | production hardening, deterministic operating procedure, release-candidate evidence | current phase |

The full milestone-by-milestone plan remains in [EXECUTION_ROADMAP_V1.md](./EXECUTION_ROADMAP_V1.md).

## Read Next

- [TECH-SPEC.md](./TECH-SPEC.md): numerical and architectural contract
- [EXECUTION_ROADMAP_V1.md](./EXECUTION_ROADMAP_V1.md): milestone sequence and validation gates
- [DETERMINISTIC_RUNBOOK.md](./DETERMINISTIC_RUNBOOK.md): reproducible build, validation, restart, and profiling workflow
- [RELEASE_CANDIDATE.md](./RELEASE_CANDIDATE.md): acceptance snapshot for the current release-candidate state
