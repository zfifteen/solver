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
-O3
-march=armv8-a
-ffast-math (optional, test impact)
```

Baseline dependencies

None beyond:

* standard library
* Apple Clang
* CMake

Profiling tools configured

* Instruments
* time-based profiling script

---

## Validation Gate

* project builds clean
* unit test framework runs
* example program executes

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

---

## Deliverables

Field allocation

Index helpers

Boundary indexing utilities

Memory alignment checks

---

## Validation Gate

Tests must verify:

* indexing correctness
* ghost cell access
* memory layout

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
(u · ∇)u
ν∇²u
```

Scheme

second-order central differencing.

---

## Deliverables

Advection kernel

Diffusion kernel

CFL diagnostic

---

## Validation Gate

Taylor-Green vortex test.

Check:

energy decay

velocity error

convergence behavior.

---

# Milestone 4 — Pressure Projection Core

## Goal

Implement fractional-step incompressible solver.

Algorithm per timestep

1. compute intermediate velocity u*
2. compute divergence of u*
3. solve pressure Poisson equation
4. correct velocity
5. enforce divergence-free field

---

## Deliverables

Projection step implementation.

Pressure RHS builder.

Velocity correction step.

---

## Validation Gate

Divergence norm must approach zero.

Test case:

static fluid

Expected result:

velocity remains zero.

---

# Milestone 5 — Linear Solver System

## Goal

Solve the pressure Poisson equation efficiently.

Solver

Conjugate Gradient

Preconditioner

Incomplete Cholesky

---

## Deliverables

Sparse matrix representation

Matrix-vector product kernel

CG iteration loop

Residual tracking

---

## Validation Gate

Poisson solve test.

Known analytic solution.

Measure:

L2 error

residual convergence rate.

---

# Milestone 6 — First Complete Simulation

## Goal

Run full incompressible solver.

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

---

## Validation Gate

Compare centerline velocity profiles against reference data.

Tolerance threshold defined.

Mass conservation must remain stable.

---

# Milestone 7 — Boundary Condition System

## Goal

Implement flexible boundary conditions.

Types

* no-slip wall
* symmetry
* inflow
* outflow
* periodic

---

## Deliverables

Boundary abstraction

BC application routines

---

## Validation Gate

Run Couette and Poiseuille flow tests.

Velocity profiles must match analytic solutions.

---

# Milestone 8 — Restart and Output

## Goal

Enable long simulation runs.

Checkpoint files contain

* velocity fields
* pressure
* timestep
* parameters

---

## Deliverables

Binary checkpoint format

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

Results must match.

---

# Milestone 9 — Verification Suite

## Goal

Prove numerical correctness.

Required tests

Grid refinement study

Time-step refinement

Mass conservation monitoring

---

## Deliverables

Automated validation scripts.

Error plots.

Convergence tables.

---

## Validation Gate

Observed order of convergence must match theoretical order.

---

# Milestone 10 — Performance Profiling

## Goal

Characterize performance on M1 Max.

Profile

pressure solver

advection kernels

memory bandwidth

thread scaling

---

## Deliverables

Performance report

cells/sec throughput

scaling curves

memory usage

---

## Validation Gate

Performance bottlenecks identified with evidence.

---

# Milestone 11 — CPU Optimization

## Goal

Optimize hotspots.

Possible improvements

loop blocking

SIMD vectorization

thread scheduling tuning

memory layout improvements

---

## Deliverables

optimized kernels

before/after benchmarks

---

## Validation Gate

measurable improvement in throughput.

---

# Milestone 12 — 3D Extension

## Goal

Extend solver to full 3D.

---

## Deliverables

3D operators

3D boundary conditions

3D Poisson solve

---

## Validation Gate

3D Taylor-Green vortex.

Correct energy decay.

---

# Milestone 13 — Metal Acceleration (Conditional)

## Only executed if profiling demands it.

Potential kernels

pressure Poisson

stencil updates

---

## Deliverables

Metal compute kernels

CPU fallback removed (M1-only still satisfied).

---

## Validation Gate

GPU implementation produces identical numerical results.

Performance gain measured.

---

# Milestone 14 — Production Hardening

## Goal

Make the solver robust.

Add

input validation

error reporting

deterministic reductions

regression tests

---

## Deliverables

CI test pipeline

regression dataset

documented benchmark suite.

---

# Milestone 15 — Release Candidate

## Final deliverables

complete documentation

verified benchmark results

performance report

reproducible test suite

restartable simulations

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
