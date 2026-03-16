# Technical Specification

**Project:** Apple-Native Incompressible Flow Solver
**Target Platform:** Apple M1 Max (ARM64)
**Scope:** Production-grade incompressible Navier–Stokes solver

---

# 1. Purpose

Develop a **numerically correct, reproducible, and maintainable solver** for the incompressible Navier–Stokes equations optimized specifically for **Apple M1 Max hardware**.

The solver must:

* produce physically correct results on standard CFD benchmarks
* run efficiently on Apple M1 Max CPU
* provide deterministic reproducible results
* support restartable long simulations
* expose profiling and diagnostics
* remain maintainable by engineers other than the author

---

# 2. Target Platform (Hard Constraint)

Supported platform:

* Hardware: **Apple M1 Max**
* CPU architecture: **ARM64**
* Operating system: **macOS (Apple Silicon)**

Unsupported:

* x86
* Linux
* Windows
* CUDA
* AMD GPUs
* cross-platform portability layers

All code and build decisions optimize specifically for **Apple Silicon memory and CPU architecture**.

---

# 3. Governing Equations

The solver implements the **incompressible Navier–Stokes equations**.

Continuity equation:

∇·u = 0

Momentum equation:

∂u/∂t + (u·∇)u = −∇p/ρ + ν∇²u + f

Where:

u : velocity field
p : pressure
ρ : density (constant)
ν : kinematic viscosity
f : body forces

---

# 4. Numerical Method

## Discretization Method

Finite Volume Method (FVM)

Reason:

* conservation properties
* robustness
* industry standard for production CFD

---

## Grid

Initial implementation:

Structured Cartesian grid

Later extension:

3D structured grids

Unstructured mesh support is **explicitly out of scope for v1**.

---

## Variable Layout

Staggered grid (MAC grid)

Reason:

* avoids pressure checkerboarding
* simpler pressure-velocity coupling
* stable discretization

---

## Spatial Discretization

Operators:

Gradient
Divergence
Laplacian

Schemes:

Advection term

* second-order central differencing initially

Diffusion term

* second-order central

Pressure gradient

* central difference

---

## Time Integration

Semi-implicit fractional step method.

Steps per timestep:

1. compute intermediate velocity u*
2. solve pressure Poisson equation
3. correct velocity
4. enforce divergence-free field

---

# 5. Pressure–Velocity Coupling

Method:

Projection method (fractional step)

Pressure Poisson equation:

∇²p = (ρ/Δt) ∇·u*

Pressure solve occurs each timestep.

---

# 6. Linear Solvers

Primary solver:

Conjugate Gradient (CG)

Preconditioner:

Incomplete Cholesky (IC)

Reason:

* robust for Poisson systems
* predictable convergence

Future extension:

Multigrid

---

# 7. Boundary Conditions

Supported boundary conditions:

Velocity Dirichlet (no-slip walls)

Velocity Neumann (symmetry)

Pressure Neumann

Pressure Dirichlet (outlet)

Periodic boundaries

Boundary conditions must be validated with dedicated tests.

---

# 8. Data Structures

Memory layout must be **cache-friendly and contiguous**.

Field storage:

velocity_x[nx][ny][nz]

velocity_y[nx][ny][nz]

velocity_z[nx][ny][nz]

pressure[nx][ny][nz]

All arrays stored as **flat contiguous buffers**.

Example:

float *velocity_x

Index mapping:

i + nx*(j + ny*k)

---

# 9. Parallelization Strategy

Primary compute model:

CPU multithreading

Threading approach:

* static domain decomposition
* thread-local work partitions
* minimal synchronization

Thread counts must be configurable.

Scaling must be tested across:

* performance cores
* efficiency cores

---

# 10. Apple-Specific Performance Strategy

Optimization targets:

Memory bandwidth utilization

Cache locality

Thread scheduling across Apple core topology

SIMD vectorization via ARM NEON where beneficial.

Profiling tools:

Instruments
Xcode performance profiler
Time-based profiling

---

# 11. GPU Acceleration (Optional Phase)

Metal acceleration considered **only if profiling proves necessary**.

Potential GPU kernels:

Pressure Poisson solver

Stencil operators

Conditions for implementation:

CPU solver fails to reach required throughput.

Metal code remains isolated in a dedicated module.

---

# 12. Input System

Configuration file:

YAML or JSON

Simulation parameters include:

grid resolution

time step

Reynolds number

boundary conditions

solver tolerances

checkpoint frequency

---

# 13. Output System

Output formats:

CSV (simple diagnostics)

Binary snapshot files (restartable)

VTK (visualization)

Output intervals configurable.

---

# 14. Restart Capability

Simulation state must support restart.

Checkpoint files contain:

velocity fields

pressure field

time step index

simulation parameters

Restart must reproduce identical continuation behavior.

---

# 15. Logging and Diagnostics

Runtime logging includes:

timestep

residuals

divergence norm

pressure solver iterations

CFL number

Logging levels:

INFO

DEBUG

ERROR

---

# 16. Validation Suite

Solver must reproduce canonical benchmarks.

Required tests:

Lid-driven cavity flow

Plane Poiseuille flow

Couette flow

Taylor–Green vortex

Manufactured solution test

Metrics evaluated:

velocity profile error

pressure field error

divergence magnitude

observed order of convergence

---

# 17. Verification Requirements

The solver must pass:

Grid refinement study

Time-step convergence

Mass conservation verification

Poisson residual convergence

Acceptance thresholds must be defined for each test.

---

# 18. Testing Infrastructure

Testing layers:

Unit tests

Numerical regression tests

Benchmark validation tests

Continuous testing must verify:

operator correctness

solver stability

output reproducibility

---

# 19. Performance Targets

Target domain sizes:

2D: up to 4096 × 4096

3D: up to 512³

Pressure solver must converge within defined tolerance.

Simulation throughput must be measured in:

cells per second

time per timestep

---

# 20. Code Organization

Proposed structure:

core
mesh representation

fields
velocity and pressure storage

operators
gradient, divergence, laplacian

solver
time stepping and projection

linsolve
linear solvers and preconditioners

bc
boundary condition implementations

io
configuration and output

tests
verification and validation

benchmarks
standard CFD cases

metal
optional GPU kernels

---

# 21. Build System

Build tool:

CMake

Compiler:

Apple Clang

Required build flags:

ARM64 optimization

vectorization enabled

debug symbol generation in development builds

---

# 22. Determinism

Solver must produce deterministic output when run with identical configuration.

Sources of nondeterminism must be controlled:

thread scheduling

floating-point reductions

random seeds

---

# 23. Documentation

Documentation must include:

equation derivations

discretization scheme

solver algorithm

benchmark results

performance profiling

known limitations

---

# 24. Known Non-Goals (v1)

The following features are **explicitly excluded**:

turbulence models

multiphase flow

compressible flow

moving meshes

adaptive mesh refinement

unstructured meshes

GPU-first execution

These may be added only after the core solver is validated.

---

# 25. Success Criteria

The solver is considered production-grade when it:

produces correct benchmark results

demonstrates numerical convergence

runs stable long simulations

restarts correctly from checkpoints

passes regression tests

operates efficiently on Apple M1 Max

is understandable and maintainable by new engineers

---
