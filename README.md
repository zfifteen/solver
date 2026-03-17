# Apple-Native Solver

Apple-Native Solver is an incompressible CFD codebase built specifically for Apple Silicon. It is designed to be numerically serious, deterministic in its reference workflow, and transparent enough to validate, restart, benchmark, and profile without turning into a portability-first research grab bag.

The project is intentionally narrow. Instead of chasing every platform and physics model at once, it focuses on a solver stack that fits Apple hardware, uses a clear numerical contract, and produces evidence you can inspect.

## Why This Project Feels Different

- Apple-native by design, not just portable to macOS after the fact
- deterministic and benchmark build profiles with an explicit reproducibility story
- structured-grid incompressible focus instead of a broad but under-validated feature surface
- validation, restart, and profiling workflows treated as first-class engineering concerns
- targeted Metal acceleration only where it is already justified by the current benchmark surface

## Supported Scope

Current supported surface:

- macOS on Apple Silicon with Apple Clang and CMake
- incompressible, constant-density flow on structured Cartesian grids
- CPU-first deterministic validation path
- lid-driven cavity, Couette, Poiseuille, and Taylor-Green verification cases
- 2D and 3D periodic Taylor-Green support
- targeted Metal backend for the 3D periodic Taylor-Green path

Intentionally out of scope for this version:

- x86, Linux, and Windows support
- unstructured meshes and adaptive mesh refinement
- compressible flow, turbulence models, and multiphase physics
- general-purpose GPU portability layers

## Proof Points

- The current validation suite reports `PASS` for operator spatial order, Taylor-Green temporal self-convergence, benchmark thresholds, and mass conservation. See [validation/latest/summary.md](validation/latest/summary.md).
- The solver has a documented deterministic run procedure plus a restartable cavity path backed by versioned checkpoints and automated continuation coverage. See [docs/DETERMINISTIC_RUNBOOK.md](docs/DETERMINISTIC_RUNBOOK.md).
- The current Apple M1 Max profiling snapshot recommends the `benchmark` profile with the default unclamped scheduler policy and reports targeted Metal speedups on the supported 3D Taylor-Green slice. See [profiling/latest/summary.md](profiling/latest/summary.md).
- The repository is in its Milestone 15 release-candidate phase for the currently supported scope. See [docs/PROJECT_STATUS.md](docs/PROJECT_STATUS.md) and [docs/RELEASE_CANDIDATE.md](docs/RELEASE_CANDIDATE.md).

## Quickstart

Configure and build the deterministic reference profile:

```bash
cmake --preset deterministic
cmake --build build/deterministic
```

Run the core automated test target:

```bash
ctest --test-dir build/deterministic --output-on-failure
```

Run one representative verification case:

```bash
build/deterministic/tools/solver_taylor_green benchmarks/taylor_green_128.cfg
```

For the fuller release-candidate workflows:

- deterministic validation suite: `./validation/run_validation_suite.py --build-dir build/deterministic --output-dir validation/latest`
- release-candidate wrapper: `tools/run_release_candidate_suite.sh`
- profiling suite: `./profiling/run_profile_suite.py --build-dir build/benchmark --output-dir profiling/latest`
- checkpoint and restart procedure: [docs/DETERMINISTIC_RUNBOOK.md](docs/DETERMINISTIC_RUNBOOK.md)

## Read Next

If you are evaluating the project:

- [docs/PROJECT_STATUS.md](docs/PROJECT_STATUS.md): current phase, implemented capabilities, limitations, and milestone snapshot
- [docs/RELEASE_CANDIDATE.md](docs/RELEASE_CANDIDATE.md): release-candidate acceptance surface and evidence links
- [validation/latest/summary.md](validation/latest/summary.md): latest validation gate results
- [profiling/latest/summary.md](profiling/latest/summary.md): latest performance snapshot

If you are implementing or extending the solver:

- [docs/TECH-SPEC.md](docs/TECH-SPEC.md): numerical and architectural contract
- [docs/EXECUTION_ROADMAP_V1.md](docs/EXECUTION_ROADMAP_V1.md): milestone-by-milestone execution plan
- [docs/DETERMINISTIC_RUNBOOK.md](docs/DETERMINISTIC_RUNBOOK.md): reproducible build, validation, restart, and profiling workflow

## Repository Shape

The codebase is organized around a few clear subsystems:

- `core/`, `operators/`, `solver/`, and `linsolve/` for the numerical implementation
- `bc/`, `io/`, and `metal/` for boundary conditions, restart/output, and the targeted GPU slice
- `tests/`, `validation/`, `profiling/`, and `tools/` for correctness, evidence generation, and measurement

## License

Licensing has not been finalized yet. Until a license file is added to the repository, do not assume one.
