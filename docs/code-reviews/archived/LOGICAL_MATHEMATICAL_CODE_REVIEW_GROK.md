# Logical and Mathematical Code Review (Grok/opencode)

**Repository:** solver
**Reviewer:** opencode (powered by grok-4.20)
**Date:** 2026-03-17
**Scope:** Logical errors, mathematical correctness in discretization, projection, operators, solvers, BCs, Metal path, validation. Focused exclusively on confirmed issues.

## Confirmed Findings

### 1. High: Drivers commit state after non-converged pressure projections (confirmed)
- **Files:**
  - `solver/lid_driven_cavity.cpp:584`
  - `solver/channel_flow.cpp:490` 
  - `solver/taylor_green.cpp:597,688`
- **Issue:** `project_velocity()` returns `ProjectionDiagnostics` with `pressure_solve.converged` flag (see `linsolve/poisson_solver.cpp:736,780`). All drivers consume only `iterations`/`relative_residual` for logging but **ignore the `converged` boolean** before updating velocity/pressure and continuing the timestep loop.
- **Mathematical impact:** Violates the incremental projection method contract. Post-projection velocity may not be divergence-free within tolerance.
- **Evidence:** Reproducible by setting `poisson_max_iterations=1` in a smoke config; run completes with high divergence_l2 (~1e-3) and residual ~0.3 while reporting success.

### 2. Medium: Validation suite L2 difference computation is unsafe
- **File:** `validation/run_validation_suite.py`
- **Functions:** `relative_l2_difference()`, `read_vtk_velocity()`, `write_taylor_temporal_reports()`
- **Issue:** `zip(left, right)` without length check. Mismatched or truncated VTK fields silently compute on prefix only (can return 0.0 erroneously).
- **Impact:** Temporal self-convergence reports and validation gates can give false positives.
- **Evidence:** `relative_l2_difference([1.0, 1000.0], [1.0]) == 0.0` (incorrect).

### 3. Medium: Metal backend precision/tolerance inconsistency
- **Files:** `metal/taylor_green_backend.mm:21-29`, `solver/taylor_green.cpp:688`
- **Issue:** `#if SOLVER_METAL_USE_FLOAT` selects float + relaxed CG tolerances (1e-4/1e-6 vs 1e-12). GPU kernels use float; final CPU `project_velocity` cleanup uses double. Config `poisson_tolerance` not consistently passed to Metal path.
- **Mathematical impact:** Cross-precision divergence in validation between CPU-only and Metal runs. Potential for non-bitwise or non-reproducible results across backends.
- **Note:** Metal path is intentionally narrow (3D periodic TG only).

## Positive Observations
- Discrete operators (`operators/discrete_operators.cpp`) implement correct 2nd-order central differences for MAC grid (gradient, divergence, laplacian).
- MGPCG (`linsolve/poisson_solver.cpp`) correctly handles null-space via mean subtraction, uses proper PCG formulas, and guards divisions.
- Ghost filling and BC mapping in `bc/` and `projection.cpp` appear consistent with projection method requirements.
- No arithmetic instabilities, incorrect stencils, or order-reduction bugs found in core math.

## Recommendations
- Add `if (!projection.pressure_solve.converged) { /* error or retry */ }` in drivers.
- Add length assertions in Python L2 helper.
- Unify tolerances between Metal and CPU paths or document the float approximation explicitly.
- Extend `tests/test_runtime.cpp` to cover convergence failure modes.

This review is documentation-only. No code changes were made.

**Previous reviews:** See `LOGICAL_MATHEMATICAL_CODE_REVIEW.md` and `CODE_REVIEW_FINDINGS.md` for additional context.