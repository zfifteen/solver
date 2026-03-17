# Execution Roadmap v1

**Project:** Apple-Native Incompressible Solver
**Platform:** Apple M1 Max (ARM64 only)

Principle: build the smallest correct system first, then extend it.

Each milestone has:

* deliverables
* validation gates
* failure signals

Progress only when the gate passes.

---

# Milestone 0 — Environment Lockdown

## Goal

Create a deterministic development environment.

## Deliverables

Repository structure

```
solver/
  core/
  operators/
  solver/
  linsolve/
  bc/
  io/
  tests/
  benchmarks/
  tools/
```

Build system

CMake configuration

Compiler flags

```
Deterministic validation profile:
  -O3
  -march=armv8-a
  no -ffast-math

Benchmark profile:
  -O3
  -march=armv8-a
  -ffast-math optional for labeled performance experiments only
```

Baseline dependencies

None beyond:

* standard library
* Apple Clang
* CMake
* Apple system frameworks as needed for threading / profiling

Profiling tools configured

* Instruments
* time-based profiling script

---

## Validation Gate

* project builds clean
* unit test framework runs
* example program executes
* deterministic and benchmark build profiles are clearly labeled

---

# Milestone 1 — Mesh and Field Infrastructure

## Goal

Establish data structures for simulation state.

## Components

Structured grid

```
struct Grid {
    nx, ny, nz
    dx, dy, dz
}
```

Field storage

```
VelocityField
PressureField
ScalarField
```

Memory layout

Flat contiguous buffers

Index mapping

```
i + nx*(j + ny*k)
```

Ghost cell support

Cell / face placement

* pressure at cell centers
* velocity components at MAC-grid faces

Precision

* `double` for all solution-state storage

Traversal contract

* `i` is unit-stride and the primary vectorization axis

---

## Deliverables

Field allocation

Index helpers

Boundary indexing utilities

Memory alignment checks

Ghost-cell fill helpers

---

## Validation Gate

Tests must verify:

* indexing correctness
* ghost cell access
* memory layout
* cell-center / face-center placement
* double-precision storage contract

---

# Milestone 2 — Discrete Operators

## Goal

Implement fundamental spatial operators.

Operators:

gradient

divergence

laplacian

All second-order finite volume stencils.

---

## Deliverables

Operator implementations

```
grad(p)
div(u)
laplacian(u)
```

Operator unit tests.

---

## Validation Gate

Manufactured solution test.

Example:

Analytic function

```
φ(x,y) = sin(x) cos(y)
```

Verify gradient and Laplacian error norms.

Expected convergence:

second-order spatial accuracy.

---

# Milestone 3 — Advection and Diffusion

## Goal

Compute RHS of momentum equation.

Term implementations

```
∇·(u u)
ν∇²u
```

Scheme

bounded second-order TVD flux-limited advection

Default limiter

van Leer

Fallback

first-order upwind for debugging only

---

## Deliverables

Advection kernel

Diffusion kernel

CFL diagnostic

Scheme-selection configuration

Limiter selection and logging

---

## Validation Gate

Taylor-Green vortex test.

Check:

energy decay

velocity error

convergence behavior.

No spurious advection oscillations in the bounded-scheme regression case.

Default validation configurations use advective CFL <= 0.5 unless a benchmark
case explicitly sets a lower value.

---

# Milestone 4 — Pressure Projection Core

## Goal

Implement second-order semi-implicit incompressible projection.

Algorithm per timestep

1. assemble explicit Adams-Bashforth advection / forcing terms
2. solve Crank-Nicolson predictor for u* via ADI Helmholtz solves
3. compute divergence of u*
4. solve pressure-correction Poisson equation with BC-consistent RHS
5. correct velocity and pressure
6. enforce divergence-free field

Startup rule

timestep 0 uses a Forward Euler bootstrap for the explicit terms before AB2 begins

---

## Deliverables

Projection step implementation.

Pressure RHS builder.

Velocity correction step.

Pressure BC mapping for each physical boundary type.

Null-space handling for pure-Neumann pressure systems.

ADI Helmholtz predictor solver with deterministic tridiagonal line solves.

---

## Validation Gate

Divergence norm must approach zero.

Test case:

static fluid

Expected result:

velocity remains zero.

Pure-Neumann cavity test must enforce zero-mean pressure and converge.

---

# Milestone 5 — Linear Solver System

## Goal

Solve the pressure Poisson equation efficiently.

Solver

MGPCG

Preconditioner

fixed V-cycle geometric multigrid with damped-Jacobi smoothing

---

## Deliverables

Matrix-free Poisson operator

Multigrid hierarchy

Documented smoother, cycle schedule, and coarse-grid solve policy

CG iteration loop

Residual tracking

---

## Validation Gate

Poisson solve test.

Known analytic solution.

Measure:

relative L2 error

relative residual convergence rate.

Acceptance:

* relative L2 error <= 1e-8
* relative residual <= 1e-10
* zero-mean pressure preserved for pure-Neumann cases

---

# Milestone 6 — First Complete Simulation

## Goal

Run full incompressible solver.

Scope note

M6 uses a hardcoded lid-driven cavity boundary-condition subset. M7 extracts this
into the general reusable boundary-condition system.

Benchmark

Lid-driven cavity flow.

Typical parameters

Re = 100

Grid

128×128

---

## Deliverables

Simulation driver

Configuration loader

Boundary condition support

Hardcoded cavity BC implementation for the benchmark case

---

## Validation Gate

Compare centerline velocity profiles against reference data.

Tolerance:

* four named centerline sample points within the accepted steady-cavity
  literature envelope, with any miss outside that envelope limited to 2
  percent

Mass conservation must remain stable with post-projection divergence L2 <= 1e-10.

---

# Milestone 7 — Boundary Condition System

## Goal

Implement flexible boundary conditions.

This milestone generalizes the cavity-specific BC subset from M6 into the full BC
system.

Types

* no-slip wall
* symmetry / slip
* prescribed inflow
* fixed-pressure outflow
* periodic

---

## Deliverables

Boundary abstraction

BC application routines

Ghost-cell fill rules

Projection-compatible pressure BC rules

---

## Validation Gate

Run Couette and Poiseuille flow tests.

Velocity profiles must match analytic solutions with relative L2 error <= 5e-3.

---

# Milestone 8 — Restart and Output

## Goal

Enable long simulation runs.

Checkpoint files contain

* velocity fields
* pressure
* timestep
* parameters
* scalar precision metadata
* format version
* checksum
* build / configuration hash

---

## Deliverables

Versioned little-endian binary checkpoint format

Checksum verification

Restart loader

VTK export

---

## Validation Gate

Restart test.

Run simulation:

```
t = 0 → 100
restart
100 → 200
```

Compare against uninterrupted run.

Results must match bitwise in the deterministic build profile.

---

# Milestone 9 — Verification Suite

## Goal

Prove numerical correctness.

Required tests

Grid refinement study

Time-step refinement

Mass conservation monitoring

Reference-data comparison for lid-driven cavity

---

## Deliverables

Automated validation scripts.

Error plots.

Convergence tables.

Named quantitative thresholds for every benchmark.

---

## Validation Gate

Observed order of convergence must be at least 1.8 where second-order accuracy is expected.

All benchmark thresholds from the technical specification must pass.

---

# Milestone 10 — Performance Profiling

## Goal

Characterize performance on M1 Max.

Profile

pressure solver

advection kernels

memory bandwidth

thread scaling

performance-core-only vs mixed P/E execution

QoS mode impact

---

## Deliverables

Performance report

cells/sec throughput

scaling curves

memory usage

recommended default thread / QoS policy

---

## Validation Gate

Performance bottlenecks identified with evidence.

Default execution mode selected from measured data, not assumption.

---

# Milestone 11 — CPU Optimization

## Goal

Optimize hotspots.

Possible improvements

loop blocking

SIMD vectorization along `i`

thread scheduling tuning

memory layout improvements

Apple Accelerate evaluation where appropriate

---

## Deliverables

optimized kernels

before/after benchmarks

---

## Validation Gate

measurable improvement in throughput.

No regression in deterministic validation results.

---

# Milestone 12 — 3D Extension

## Goal

Extend solver to full 3D.

---

## Deliverables

3D operators

3D boundary conditions

3D Poisson solve

validated 256^3 execution path

---

## Validation Gate

3D Taylor-Green vortex.

Correct energy decay.

512^3 remains a stretch target until profiling and optimization gates pass.

---

# Milestone 13 — Metal Acceleration (Conditional)

## Executed once profiling showed the 3D Taylor-Green path was the next meaningful acceleration target.

Current scoped backend:

3D periodic Taylor-Green only

explicit backend switch on the existing Taylor-Green executable

CPU path retained as the reference path for comparison and fallback

---

## Deliverables

targeted Metal compute backend in a dedicated module

CPU fallback retained for validation and reproducibility comparison

CPU-vs-Metal comparison report in the profiling harness

---

## Validation Gate

Metal path passes the supported 3D Taylor-Green benchmark gate.

CPU-vs-Metal velocity / energy agreement measured on a small 3D case.

Performance gain measured on the 3D `256^3` smoke path.

---

# Milestone 14 — Production Hardening

## Goal

Make the solver robust.

Add

input validation

error reporting

deterministic reductions

regression tests

checkpoint compatibility tests

manual deterministic benchmark-threshold gate

---

## Deliverables

fast CI test pipeline

regression dataset

documented benchmark suite and manual threshold gate.

---

# Milestone 15 — Release Candidate

## Final deliverables

complete documentation

verified benchmark results

performance report

reproducible test suite

restartable simulations

documented deterministic build and run procedure

---

# Expected Development Sequence

Approximate order:

```
M0  Environment
M1  Grid & fields
M2  Operators
M3  Advection/diffusion
M4  Projection step
M5  Linear solver
M6  First benchmark
M7  Boundary conditions
M8  Restart
M9  Verification
M10 Profiling
M11 Optimization
M12 3D
M13 Metal acceleration
M14 Hardening
M15 Release
```

---

# Development Rule

At every stage:

If validation fails, stop and fix before adding new features.

Production solvers are built by **eliminating uncertainty step by step**, not by adding complexity quickly.

---
