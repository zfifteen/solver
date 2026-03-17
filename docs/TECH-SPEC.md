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

∂u/∂t + ∇·(u⊗u) = −∇p/ρ + ν∇²u + f

Where:

u : velocity field
p : pressure
ρ : density (constant)
ν : kinematic viscosity
f : body forces

For constant density incompressible flow, the conservative convective form above is
analytically equivalent to (u·∇)u when ∇·u = 0. The v1 finite-volume
discretization uses the conservative form so that face-flux balances remain
discretely conservative before projection.

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
Convective flux divergence

Schemes:

Advection term

* conservative divergence form ∇·(u⊗u)
* bounded second-order TVD flux-limited face-flux scheme
* default limiter: van Leer
* limiter choice must be configurable and logged
* first-order upwind mode permitted only as a debugging fallback
* pure central differencing is not an accepted production default

Diffusion term

* second-order central

Pressure gradient

* central difference

Interpolation between cell centers and faces must be second-order consistent with
the MAC grid arrangement.

---

## Time Integration

Second-order semi-implicit incremental pressure-correction method.

Term treatment:

* advection: explicit Adams-Bashforth 2
* diffusion: implicit Crank-Nicolson
* pressure: projection every timestep
* timestep control: fixed Δt in v1; adaptive time stepping is out of scope
* AB2 startup: timestep 0 uses a Forward Euler bootstrap for the explicit terms,
  then switches to AB2 from timestep 1 onward
* Crank-Nicolson diffusion solve: ADI factorization of the structured-grid
  Helmholtz operator with deterministic tridiagonal line solves

Timestep selection must satisfy a configured advective CFL limit and the selected
Δt must be logged for every run.

Validation-harness default:

* advective CFL <= 0.5 unless a benchmark definition explicitly requires a lower
  value

Steps per timestep:

1. assemble explicit convective and body-force terms
2. solve the implicit predictor for intermediate velocity u*
3. solve the pressure-correction Poisson equation for phi with BC-consistent RHS
4. correct velocity and update the total pressure
5. enforce divergence-free velocity and log post-projection norms

---

# 5. Pressure–Velocity Coupling

Method:

Incremental projection method (fractional step)

Pressure-correction Poisson equation:

∇²phi = (ρ/Δt) ∇·u*

Pressure solve occurs each timestep.

Pressure boundary-condition rules:

* the projection solve is for the pressure-correction variable phi, not the
  total physical pressure field
* no-slip walls, symmetry planes, and prescribed-velocity inflows use
  homogeneous Neumann correction BCs
* prescribed-pressure outlets use homogeneous Dirichlet correction BCs
* periodic boundaries remain periodic for both pressure correction and velocity
* closed domains with only Neumann correction boundaries enforce solvability by
  removing the mean of the Poisson RHS and solving with a zero-mean correction
  constraint

The projection implementation must apply boundary conditions and null-space
treatment consistently in the predictor, Poisson RHS assembly, and velocity
correction steps.

---

# 6. Linear Solvers

Primary solver:

Geometric multigrid-preconditioned Conjugate Gradient (MGPCG)

Preconditioner:

Fixed V-cycle geometric multigrid with damped-Jacobi smoothing

Reason:

* scalable for large structured Poisson systems
* compatible with CPU multithreading
* required to support the 3D target domain sizes

Implementation requirements:

* matrix-free structured-grid stencil operator
* damped-Jacobi smoother as the deterministic default
* fixed V-cycle schedule
* deterministic residual and dot-product reductions
* direct coarse-grid solve on the coarsest level
* reference plain CG permitted only for small test problems

Predictor diffusion solver requirements:

* the Crank-Nicolson predictor Helmholtz systems use ADI splitting on the
  structured grid
* each ADI substep uses deterministic tridiagonal line solves
* predictor-solve tolerances and iteration structure must be identical across
  deterministic reruns

Incomplete Cholesky is out of scope for v1 because it conflicts with the parallel
execution strategy.

---

# 7. Boundary Conditions

Supported boundary conditions:

Velocity Dirichlet (no-slip walls, prescribed inflow)

Velocity symmetry / slip

Pressure Neumann

Pressure Dirichlet (outlet)

Periodic boundaries

Ghost-cell convention:

* one ghost-cell layer minimum around every field
* pressure stored at cell centers
* velocity components stored on faces of the MAC grid
* wall-normal velocity imposed directly at boundary faces
* tangential wall velocities filled through ghost values consistent with no-slip
  or symmetry constraints
* pressure ghost values filled from the prescribed Dirichlet value or the
  BC-specific normal derivative

Boundary conditions must be defined as a complete mapping from physical BC to:

* velocity value/gradient treatment
* pressure value/gradient treatment
* ghost-cell fill rule
* projection correction rule

All boundary conditions must be validated with dedicated tests.

---

# 8. Data Structures

Memory layout must be **cache-friendly, contiguous, and deterministic**.

Field storage:

velocity_x[nx][ny][nz]

velocity_y[nx][ny][nz]

velocity_z[nx][ny][nz]

pressure[nx][ny][nz]

All arrays stored as **flat contiguous buffers**.

Precision policy:

* all solution state is stored in `double`
* pressure-solver residuals, dot products, and norms are accumulated in `double`
* mixed-precision production paths are out of scope for v1 reference results

Example:

double *velocity_x

Index mapping:

i + nx*(j + ny*k)

Layout requirements:

* the `i` dimension is unit-stride and is the primary SIMD/vectorization axis
* loop nests must traverse `i` in the innermost position unless profiling proves
  another order is superior
* ghost-cell storage is part of each field allocation, not a side structure

---

# 9. Parallelization Strategy

Primary compute model:

CPU multithreading

Threading approach:

* fixed-size worker pool
* deterministic domain decomposition with fixed reduction trees
* performance-core-first execution by default
* mixed performance/efficiency-core execution allowed only with explicit weighted
  partitions and profiling evidence
* minimal synchronization at timestep phase boundaries

Thread counts must be configurable.

Thread QoS class must be configurable and recorded in logs and checkpoints.

Scaling must be tested across:

* performance cores
* efficiency cores

---

# 10. Apple-Specific Performance Strategy

Optimization targets:

Memory bandwidth utilization

Cache locality with unit-stride `i` sweeps

Thread scheduling across Apple core topology

SIMD vectorization via ARM NEON where beneficial.

Platform-specific execution rules:

* compute threads must run at an explicit QoS class
* the default benchmark and production profile uses `QOS_CLASS_USER_INITIATED`
* reports must state whether a run used performance-core-only or mixed P/E mode
* Apple Accelerate / vecLib may be used for BLAS-like kernels if validation and
  deterministic-build requirements remain satisfied

Profiling tools:

Instruments
Xcode performance profiler
Time-based profiling

---

# 11. GPU Acceleration (Optional Phase)

Metal acceleration is opened only after profiling establishes that the 3D CPU path is the next meaningful optimization target.

Current v1 GPU scope:

* Metal backend is isolated in a dedicated `metal/` module
* supported GPU surface is intentionally narrow: 3D periodic Taylor-Green only
* backend selection is explicit through Taylor-Green config / CLI backend selection
* CPU path remains the validation source of truth for the general solver

Current implementation notes:

* Apple Metal on the target machine does not support `double` kernels, so the first GPU slice uses float working storage on device
* final result publication performs a CPU projection cleanup before finalization so the reported state still satisfies the strict divergence contract
* cavity, channel, restart, and general-BC paths remain CPU-only in this milestone

---

# 12. Input System

Configuration file:

YAML or JSON

Simulation parameters include:

grid resolution

time step

time integration scheme

advection scheme / limiter mode

Reynolds number

boundary conditions

solver tolerances

thread count and QoS mode

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

scalar precision metadata

format version

integrity checksum

build / configuration hash

Binary checkpoint requirements:

* explicit versioned format
* defined little-endian encoding
* checksum verification before restart
* forward-compatible handling of added metadata fields

Restart must reproduce identical continuation behavior in the deterministic build
profile when executable, thread configuration, and input configuration match.

---

# 15. Logging and Diagnostics

Runtime logging includes:

timestep

residuals

divergence norm

pressure solver iterations

multigrid cycles

pressure null-space correction status

CFL number

active advection scheme

active limiter

thread count and QoS mode

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

kinetic energy error for unsteady benchmarks

All acceptance thresholds in section 17 are normative.

---

# 17. Verification Requirements

The solver must pass:

* Manufactured solution: observed spatial order >= 1.8 for gradient,
  divergence, and Laplacian under uniform 2x mesh refinement.
* Time-step convergence: observed temporal order >= 1.8 on the selected
  second-order time integrator.
* Couette flow: relative L2 velocity-profile error <= 5e-3 on the reference
  128 x 128 case.
* Plane Poiseuille flow: relative L2 velocity-profile error <= 5e-3 on the
  reference 128 x 128 case.
* Lid-driven cavity, Re = 100, 128 x 128: four named centerline sample points
  in the validation harness must lie within the accepted steady-cavity
  literature envelope, with any miss outside that envelope limited to 2
  percent.
* Taylor-Green vortex: normalized kinetic-energy error <= 1 percent over the
  benchmark time horizon.
* Poisson analytic solve: relative L2 pressure error <= 1e-8 and relative
  residual <= 1e-10 in double precision.
* Mass conservation: post-projection divergence L2 norm <= 1e-10 on 2D
  reference cases.
* Restart reproducibility: restarted and uninterrupted deterministic-profile
  runs must match bitwise.

Each benchmark case must name its reference dataset, norm definition, sample
times, and grid / timestep parameters in the validation harness.

Unless a benchmark definition explicitly overrides it, the validation harness must
use the default advective CFL ceiling defined in section 4.

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

2D validated target: up to 4096 × 4096

3D validated target: up to 256³

3D stretch target: 512³ after multigrid, profiling, and optimization gates pass

Pressure solver must converge to a relative residual of at least 1e-10, or to a
stricter case-specific tolerance if required by validation.

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

Build profiles:

* deterministic validation profile: no `-ffast-math`; stable floating-point
  semantics required
* benchmark profile: `-O3`; `-ffast-math` allowed only for clearly labeled
  performance experiments and never for reference-result generation

---

# 22. Determinism

Solver must produce deterministic output when run with identical configuration.

Sources of nondeterminism must be controlled:

* thread count
* QoS class
* domain partitioning
* floating-point reductions
* random seeds
* restart metadata / build configuration

Deterministic execution requirements:

* fixed reduction order for dot products, norms, and global diagnostics
* no schedule-dependent atomics in the numerical path
* deterministic build profile is the source of truth for validation and restart
  equivalence
* any non-deterministic benchmark profile must be clearly labeled as such

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

adaptive time stepping

unstructured meshes

GPU-first execution

mixed-precision production solves

These may be added only after the core solver is validated.

---

# 25. Success Criteria

The solver is considered production-grade when it:

meets all quantitative validation thresholds

demonstrates the expected spatial and temporal convergence

runs stable long simulations

restarts correctly from versioned checkpoints

passes deterministic regression tests

operates efficiently on Apple M1 Max

documents benchmark references, build modes, and known limitations

is understandable and maintainable by new engineers

---
