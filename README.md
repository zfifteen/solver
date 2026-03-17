# Apple-Native Solver

This repository is the beginning of a very specific kind of CFD project: an incompressible flow solver that is meant to feel native to Apple Silicon rather than merely portable to it. The long-term goal is not to assemble a generic research code that happens to compile on a Mac. The goal is to build a solver that is numerically serious, reproducible, restartable, profileable, and understandable, while also being unapologetically tuned for the realities of Apple hardware and toolchains.

That combination matters here. A lot of scientific software grows by accumulating features first and worrying about rigor later. This project is trying to do the opposite. The solver is being designed from the outset around a clear numerical method, explicit validation gates, deterministic execution rules, and a roadmap that only moves forward when the previous layer is correct. If the end result works, it should be the kind of codebase a new engineer can open, reason about, benchmark, and trust.

Right now, the repository is still early, but it has moved well past the "just scaffolding" stage. It now contains multiple full simulation and verification paths rather than only isolated numerical building blocks. The codebase has the locked-down build environment, the structured-grid and field-storage layer, the discrete operators, the transport-term kernels, the projection machinery, the matrix-free MGPCG pressure solver, the cavity benchmark path, the generalized boundary-condition system exercised by analytic Couette and Poiseuille verification cases, the restart/output path for the first full driver, the automated broader verification suite with convergence reports, the profiling layer with microbenchmarks and Instruments-backed hotspot summaries, and now a first targeted Metal backend for the 3D periodic Taylor-Green path. In other words, this repo is already opinionated about how the solver should be built, validated, measured, and evolved before the deeper production-hardening passes arrive.

If you are reading this as a developer, the shortest useful summary is: this project is building toward a production-grade, Apple-Silicon-native incompressible flow solver, and the repository currently reflects Milestone 13 of that plan.

## Current Status

The repository is currently at **Milestone 13: Targeted Metal Acceleration**.

Implemented today:

- Apple-Silicon-only CMake configuration using Apple Clang
- two labeled build profiles: `deterministic` and `benchmark`
- validated `Grid` representation for structured Cartesian domains
- MAC-aware pressure, scalar, and velocity field storage types
- flat contiguous `double` buffers with explicit ghost-cell-aware indexing
- boundary slab and ghost-layer helpers for future BC work
- second-order discrete `gradient`, `divergence`, and `laplacian` operators
- manufactured-solution convergence checks for the operator layer
- transport-term kernels for advection and diffusion
- CFL diagnostics plus explicit advection-scheme and limiter configuration
- bounded TVD advection with `van Leer`, a first-order upwind fallback, and an optional central scheme for controlled investigations
- pressure-projection path with BC mapping and null-space handling
- deterministic ADI predictor solve with tridiagonal line solves
- matrix-free pressure Poisson operator
- MGPCG pressure solve with fixed V-cycle geometric multigrid preconditioning
- damped-Jacobi multigrid smoothing and a direct coarse-grid solve
- residual-history tracking plus multigrid policy diagnostics
- full hardcoded lid-driven cavity timestepper for the Milestone 6 path
- generalized physical-boundary abstraction in `bc/boundary_conditions.hpp`
- shared velocity, pressure-correction, and total-pressure ghost-fill routines
- deterministic cavity benchmark executable plus smoke and validation configs
- literature-envelope validation for the Re = 100 cavity gate
- channel-flow verification driver covering Couette and Poiseuille cases
- deterministic channel benchmark executable plus smoke and 128 x 128 validation configs
- analytic profile validation for Couette and Poiseuille at the M7 acceptance threshold
- stateful cavity stepping API for deterministic continuation and restart tests
- versioned little-endian binary checkpoint format with checksum verification
- build/configuration hash validation on checkpoint load
- checkpoint-based cavity restart loader with bitwise deterministic continuation coverage
- legacy VTK export for cell-centered pressure and reconstructed velocity
- Taylor-Green vortex verification driver with analytic kinetic-energy decay validation
- dimension-aware Taylor-Green support covering both 2D and 3D periodic cases
- manufactured-solution operator-verification tool for reusable spatial convergence studies
- automated validation harness that writes convergence tables, benchmark summaries, mass-conservation summaries, and SVG plots under `validation/latest/`
- generated Milestone 9 evidence artifacts for operator spatial convergence, Taylor-Green temporal self-convergence, and deterministic benchmark thresholds
- dedicated advection and pressure-solve microbenchmark executables
- automated profiling harness that writes kernel, throughput, policy-study, thread-scaling, and hotspot reports under `profiling/latest/`
- Instruments-backed hotspot summaries plus sampled P-core / E-core share reporting on Apple Silicon
- measured default scheduler-policy recommendation for the current M1 Max benchmark path
- optimized unchecked storage-access path for the hottest structured-field kernels
- measured Milestone 10 to Milestone 11 before/after benchmark comparison under `profiling/latest/`
- 3D Taylor-Green smoke and validation configs plus a short `256^3` execution-path smoke config
- targeted Metal backend dispatch for the 3D periodic Taylor-Green executable
- dedicated `metal/` module with Objective-C++ host orchestration and compiled Metal shaders
- explicit `backend=cpu|metal` Taylor-Green configuration and `--backend` CLI override
- targeted Metal comparison reporting in the profiling harness
- a smoke-test executable that reports build/runtime metadata
- a minimal test executable wired into CTest
- a simple time-based profiling helper script
- scaffolded module layout for later solver milestones beyond the infrastructure layer

What is not implemented yet:

- broader 3D full-case coverage beyond the Taylor-Green path
- deeper solver optimization and true multithreaded scaling beyond the current first M11 pass
- broader 3D and higher-Reynolds-number verification beyond the current reference suite

Current implementation note:

- the current cavity and channel drivers are still intentionally 2D-only even though the Taylor-Green path is now 3D-capable
- the current full-case drivers remain intentionally specialized even though they now share the generalized boundary-condition system underneath
- the cavity benchmark validates against a steady-cavity literature envelope rather than treating the coarse Ghia point table as the sole source of truth
- the current restart/output path is implemented for the cavity driver first, which is the milestone-appropriate full simulation path
- the current broader verification suite is deterministic-profile-first and intentionally centered on the milestone reference cases rather than a large benchmark catalog
- the current profiling recommendation is based on measured Apple M1 Max runs and currently favors the `benchmark` build with the default unclamped scheduler policy
- the current M11 pass improves kernel and end-to-end throughput substantially without changing the numerical method
- the current solver path is still effectively single-threaded, so this milestone closes CPU optimization before deeper thread-scaling work
- the current Taylor-Green Metal path is intentionally narrow: 3D, fully periodic, Taylor-Green-only
- the current Metal backend uses float working storage on GPU because Apple Metal does not support `double` kernels on this machine, and it performs a final CPU projection cleanup before result finalization to restore the strict divergence contract
- the current cavity and channel drivers remain CPU-only

## Project Goals

The project is working toward a solver that can:

- solve the incompressible Navier-Stokes equations on structured Cartesian grids
- produce benchmark-correct results on standard CFD validation cases
- run efficiently on Apple M1 Max class hardware
- produce deterministic reference results under a locked validation profile
- support restartable long-running simulations
- expose diagnostics and profiling information instead of treating performance as an afterthought
- remain maintainable for engineers who did not originally author it

## Platform and Scope

This project is intentionally narrow in scope.

- Supported platform: macOS on Apple Silicon, with Apple Clang and CMake
- Primary hardware target: Apple M1 Max
- Primary execution path: CPU-first, multithreaded, Apple-native
- Initial mesh scope: structured Cartesian grids
- Initial physics scope: incompressible, constant-density flow

Explicitly out of scope for the v1 solver:

- x86, Linux, and Windows support
- CUDA, AMD GPU paths, and portability layers
- unstructured meshes
- adaptive mesh refinement
- compressible flow
- turbulence models
- multiphase flow

This narrow scope is deliberate. The project is optimizing for correctness and a clean implementation path before it optimizes for breadth.

## Numerical Direction

The full solver architecture is described in [TECH-SPEC.md](TECH-SPEC.md), but the current direction is already fixed:

- discretization: finite volume method on a staggered MAC grid
- convective form: conservative `∇·(u⊗u)` form
- default advection scheme: bounded second-order TVD with the `van Leer` limiter
- optional advection alternatives: first-order upwind and centered fluxes for controlled comparisons
- time integration: second-order semi-implicit pressure-correction scheme
- explicit term startup: Forward Euler bootstrap into Adams-Bashforth 2
- diffusion treatment: Crank-Nicolson with ADI Helmholtz line solves
- pressure solve: MGPCG with fixed V-cycle multigrid and damped-Jacobi smoothing
- precision policy: `double` for solution state and pressure-solver reductions
- validation default: advective CFL `<= 0.5` unless a benchmark case says otherwise
- current full-case benchmarks: lid-driven cavity at `Re = 100`, analytic Couette and Poiseuille profile preservation cases, and both 2D and 3D Taylor-Green vortex decay

Parts of that numerical path are now implemented at the operator, transport, projection, linear-solver, boundary-condition, and restart/output layers, and the rest remains the governing design contract for the upcoming milestones.

## Repository Layout

The repository is organized around the milestone structure defined in the roadmap:

```text
solver/
  core/
  operators/
  solver/
  linsolve/
  metal/
  bc/
  io/
  tests/
  benchmarks/
  tools/
```

What those directories mean in practice right now:

- `core/`: runtime/build metadata plus the first structured-grid and field-storage layer
- `operators/`: second-order structured-grid discrete operators
- `linsolve/`: matrix-free pressure operator and MGPCG linear-solver implementation
- `solver/`: transport-term kernels, projection helpers, the cavity driver, the channel-flow verification driver, and the Taylor-Green verification driver
- `metal/`: the targeted Metal backend for the 3D periodic Taylor-Green path
- `tools/`: the smoke executable, cavity benchmark executable, channel verification executable, Taylor-Green executable, operator-verification executable, kernel microbenchmark executables, and profiling helper
- `tests/`: the CTest-backed validation executable
- `benchmarks/`: cavity plus channel-flow smoke and validation configs
- `bc/`: the shared boundary-condition abstraction used by the solver layer
- `io/`: checkpoint serialization, restart loading, and VTK export
- `validation/`: the automated Milestone 9 harness plus the latest generated reports and plots
- `profiling/`: the automated Milestone 11 harness, stored M10 baseline data, and the latest generated profiling reports and plots

The technical and roadmap documents live at the repository root and currently act as the primary design references.

## Build Profiles

The build is intentionally split into two named profiles:

- `deterministic`: the reference build for validation and reproducibility work
- `benchmark`: the performance-oriented build for clearly labeled experiments

At the current stage, those profiles differ mainly in floating-point behavior:

- `deterministic` disables `-ffast-math`
- `benchmark` enables `-ffast-math` when the compiler accepts it

Both profiles target Apple Silicon and require Apple Clang.

## Building

Configure and build the deterministic profile:

```bash
cmake --preset deterministic
cmake --build build/deterministic
```

Configure and build the benchmark profile:

```bash
cmake --preset benchmark
cmake --build build/benchmark
```

The current build produces:

- `solver_core`: a small core library for runtime/build metadata and mesh/field infrastructure
- `solver_operators`: the discrete-operator library built on top of the core field layer
- `solver_linsolve`: matrix-free Poisson application plus MGPCG / geometric multigrid pressure solve
- `solver_momentum`: advection, diffusion, CFL, projection utilities, the cavity driver, and the channel-flow verification driver layered on top of the linsolve path
- `solver_io`: checkpoint/restart and VTK export layered on top of the current solver state types
- `solver_metal`: the targeted Metal backend library for 3D periodic Taylor-Green runs
- `solver_example`: a smoke-test executable under `build/<profile>/tools/`
- `solver_cavity`: the lid-driven cavity benchmark executable under `build/<profile>/tools/`
- `solver_channel`: the Couette / Poiseuille verification executable under `build/<profile>/tools/`
- `solver_taylor_green`: the Taylor-Green vortex verification executable under `build/<profile>/tools/`
- `solver_operator_mms`: the manufactured-solution operator verification executable under `build/<profile>/tools/`
- `solver_advection_profile`: the advection microbenchmark executable under `build/<profile>/tools/`
- `solver_pressure_profile`: the pressure-solver microbenchmark executable under `build/<profile>/tools/`
- `solver_tests`: a minimal test executable under `build/<profile>/tests/`

## Testing

Run the deterministic-profile tests:

```bash
ctest --test-dir build/deterministic --output-on-failure
```

Run the benchmark-profile tests:

```bash
ctest --test-dir build/benchmark --output-on-failure
```

Run the cavity smoke config:

```bash
build/deterministic/tools/solver_cavity benchmarks/lid_driven_cavity_smoke.cfg
```

Run the deterministic Milestone 6 benchmark gate:

```bash
build/deterministic/tools/solver_cavity benchmarks/lid_driven_cavity_re100_128.cfg
```

Run the deterministic Milestone 7 benchmark gates:

```bash
build/deterministic/tools/solver_channel benchmarks/channel_couette_128.cfg
build/deterministic/tools/solver_channel benchmarks/channel_poiseuille_128.cfg
```

Run the deterministic 2D Taylor-Green benchmark gate:

```bash
build/deterministic/tools/solver_taylor_green benchmarks/taylor_green_128.cfg
```

Run the deterministic Milestone 12 3D Taylor-Green benchmark gate:

```bash
build/deterministic/tools/solver_taylor_green benchmarks/taylor_green_3d_64.cfg
```

Run the targeted Metal version of that same 3D gate:

```bash
build/deterministic/tools/solver_taylor_green benchmarks/taylor_green_3d_64.cfg --backend metal
```

Run the short `256^3` execution-path smoke:

```bash
build/deterministic/tools/solver_taylor_green benchmarks/taylor_green_3d_256_smoke.cfg
```

Run the short `256^3` targeted Metal smoke:

```bash
build/benchmark/tools/solver_taylor_green benchmarks/taylor_green_3d_256_smoke.cfg --backend metal
```

Generate the full deterministic validation report set:

```bash
./validation/run_validation_suite.py --build-dir build/deterministic --output-dir validation/latest
```

Generate the full Milestone 11 profiling and comparison report set:

```bash
./profiling/run_profile_suite.py --build-dir build/benchmark --output-dir profiling/latest
```

Write a checkpoint and VTK file from a short cavity run:

```bash
build/deterministic/tools/solver_cavity benchmarks/lid_driven_cavity_smoke.cfg --steps 6 --checkpoint-out /tmp/solver.chk --vtk-out /tmp/solver.vtk
```

Restart from that checkpoint for more steps:

```bash
build/deterministic/tools/solver_cavity benchmarks/lid_driven_cavity_smoke.cfg --restart /tmp/solver.chk --steps 6
```

Today’s GitHub Actions coverage is intentionally limited to the fast deterministic build-and-test gate. The current test executable checks that:

- the build profile is one of the locked profile names
- the runtime platform is Apple Silicon
- the generated build banner includes the active profile
- grid dimensions, spacings, and coordinate helpers behave as expected
- pressure-field indexing matches the flat-buffer mapping contract
- ghost layers and boundary slab helpers address the expected storage regions
- storage is aligned and unit-stride in the `i` direction
- MAC-grid cell-center and face-center placement is correct
- field storage uses `double`
- manufactured-solution error norms for gradient, divergence, and Laplacian decrease under refinement
- observed operator convergence is at least second-order in the Milestone 2 validation case
- diffusion matches the expected scaled Laplacian behavior on a smooth field
- Taylor-Green step behavior improves under refinement
- the bounded TVD advection regression case avoids overshoot and undershoot
- CFL diagnostics and advection-config labels report the expected values
- physical BC types map to the expected projection-side pressure BC rules
- the ADI predictor preserves the quiescent state
- the periodic ADI sweep preserves an exact Couette state
- the ADI predictor matches a dense factorized reference solve on a nontrivial case
- the generalized total-pressure ghost fill follows the normal-momentum balance and fixed-pressure contracts
- the cavity predictor RHS includes the factorization correction term
- the static-fluid projection path keeps the zero solution unchanged
- the pure-Neumann projection path removes the pressure null space and restores a divergence-free field
- the MGPCG solve recovers a known discrete Dirichlet solution to the required relative error
- the MGPCG solve preserves zero-mean pressure for a pure-Neumann problem
- the MGPCG residual converges to the required relative tolerance and reports the fixed policy metadata
- the Couette and Poiseuille verification paths load configs, apply the intended BC sets, and pass their analytic profile gates
- the Taylor-Green verification path loads both 2D and 3D configs, applies the intended periodic/symmetry BC sets, and passes its decay gates
- the Taylor-Green backend parser accepts explicit CPU vs Metal selection and rejects unsupported Metal cases
- the targeted Metal Taylor-Green path stays close to the CPU reference on a small 3D comparison case while preserving the analytic benchmark gate on the supported 3D benchmark surface
- checkpoint round-trips preserve the cavity continuation state bitwise
- corrupted checkpoints fail checksum verification before restart
- split cavity runs reproduce uninterrupted deterministic runs bitwise after restart
- VTK export writes the expected legacy dataset structure and field names
- the cavity smoke run remains stable and divergence controlled
- the cavity validation harness samples the named reference points correctly and rejects out-of-range results

The broader verification harness that closes Milestone 9 lives outside `solver_tests` in `validation/run_validation_suite.py`. It is run manually when we want full numerical evidence rather than on every GitHub Actions run. Its generated outputs in `validation/latest/` currently include:

- operator spatial-convergence tables and SVG plots
- Taylor-Green temporal self-convergence tables and SVG plots
- deterministic benchmark-threshold summaries for Couette, Poiseuille, lid-driven cavity, and both 2D and 3D Taylor-Green cases
- mass-conservation summaries across the reference cases

Milestone 11 profiling and optimization evidence also lives outside `solver_tests`, in `profiling/run_profile_suite.py`. Its generated outputs in `profiling/latest/` currently include:

- hardware and platform summaries for the active Apple Silicon machine
- advection and pressure-solver kernel microbenchmark tables and throughput plots
- end-to-end throughput summaries for the current deterministic reference cases
- QoS / scheduler-policy studies with sampled P-core and E-core shares
- Time Profiler hotspot summaries grouped by solver subsystem
- explicit comparison tables against the stored Milestone 10 M1 Max baseline
- a baseline thread-scaling report that makes the current single-threaded state explicit

## Profiling

The repository now has two profiling layers:

```bash
tools/profile_example.sh build/deterministic
```

That helper still gives a quick `/usr/bin/time -lp` baseline for the smoke executable.

For the real Milestone 11 profiling and comparison pass, use:

```bash
./profiling/run_profile_suite.py --build-dir build/benchmark --output-dir profiling/latest
```

That harness runs kernel microbenchmarks, end-to-end benchmark timing, QoS policy studies, Instruments Time Profiler captures, and a direct comparison against the stored Milestone 10 M1 Max baseline. On the current Apple M1 Max reference machine, the measured recommendation is still to use the `benchmark` profile with the default unclamped scheduler policy for performance experiments. The current compute-thread recommendation remains `1`, because this first M11 pass focused on hot-kernel efficiency before deeper thread-scaling work.

If you are running that long profiling pass from an attached terminal UI and want to avoid keeping the session live for the whole `xctrace` run, use:

```bash
tools/profile_suite_background.sh build/benchmark profiling/latest
```

That helper detaches the profiling suite, writes a timestamped log under `profiling/latest/`, and records the active PID in `profiling/latest/run_profile_suite.pid`.

## Roadmap

The implementation plan is spelled out in [EXECUTION_ROADMAP_V1.md](EXECUTION_ROADMAP_V1.md). The short version looks like this:

1. Milestone 0: lock the environment, build system, and validation scaffolding
2. Milestone 1: add grid, field, indexing, and ghost-cell infrastructure
3. Milestone 2: implement core discrete operators
4. Milestone 3: implement advection and diffusion terms
5. Milestone 4: implement the projection step
6. Milestone 5: implement the pressure linear solver system
7. Milestone 6: hardcoded cavity benchmark validation
8. Milestone 7: boundary-condition generalization plus Couette and Poiseuille verification
9. Milestone 8: restart/output for the current full simulation path
10. Milestone 9: broader verification with convergence studies and benchmark reporting
11. Milestone 10: performance profiling and baseline hotspot analysis
12. Milestone 11: CPU optimization with before/after benchmark comparison
13. Milestone 12: 3D Taylor-Green validation plus a short `256^3` execution-path proof
14. Milestone 14 and beyond: production hardening, a fast CI gate plus a manual deterministic threshold gate, broader scaling work, and expansion beyond the first targeted GPU slice

The repository has completed Milestone 13 with a deliberately narrow Metal backend for 3D periodic Taylor-Green. The next work is no longer deciding whether to open Metal at all; it is deciding how far to generalize the GPU path without diluting the solver’s validation discipline.

For Milestone 14, GitHub Actions is intentionally limited to fast deterministic build-and-test checks. The full deterministic benchmark-threshold gate is a manual developer workflow run through `validation/run_validation_suite.py`, not a default CI job.

The important project rule is simple: **do not advance to the next milestone unless the current validation gate passes**.

## Reference Documents

The two main design documents in this repository are:

- [TECH-SPEC.md](TECH-SPEC.md): the numerical, architectural, determinism, and validation contract for the solver
- [EXECUTION_ROADMAP_V1.md](EXECUTION_ROADMAP_V1.md): the milestone-by-milestone implementation sequence and validation gates

If you need to understand why the repository is structured the way it is, start with those two files.

## License

Licensing has not been finalized yet. Until a license file is added to the repository, do not assume one.
