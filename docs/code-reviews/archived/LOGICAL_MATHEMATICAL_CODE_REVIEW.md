# Logical and Mathematical Code Review

**Repository:** `solver`  
**Reviewer:** Codex  
**Date:** 2026-03-17

## Scope and Method

This review focused only on confirmed logical or mathematical errors in the current `main`-derived codebase. I inspected the pressure solver, projection path, advection discretization, MAC-grid operators, driver update loops, and validation harness. I also ran `ctest --preset deterministic`, which passed, and used one targeted runtime repro to confirm a driver-level failure mode.

## Confirmed Findings

### 1. High: CPU driver loops commit states after non-converged pressure projections

- **Affected files / functions:**
  - `solver/channel_flow.cpp` - channel-flow timestep loop
  - `solver/lid_driven_cavity.cpp` - `advance_lid_driven_cavity_impl(...)`
  - `solver/taylor_green.cpp` - `run_taylor_green_steps(...)`
  - `solver/taylor_green.cpp` - final CPU cleanup path inside `run_taylor_green(...)` for Metal runs
- **Why this is wrong:**
  - `project_velocity(...)` returns `ProjectionDiagnostics`, including `pressure_solve.converged`, `pressure_solve.iterations`, and `pressure_solve.relative_residual`.
  - The callers use the returned `corrected` velocity and update `pressure_total` immediately, but none of these paths check `pressure_solve.converged` before committing the new state.
  - That means the solver can advance with a projection that failed its own tolerance contract, which is mathematically inconsistent with the incompressibility update the drivers intend to apply.
- **Impact:**
  - A run can continue and publish apparently valid metrics even when the pressure solve did not reach the requested tolerance.
  - This can leave the state materially under-projected and inflate post-step divergence while still looking like a normal completed run.
- **Concrete evidence:**
  - Code inspection shows the returned projection diagnostics are used for reporting only; the convergence bit is ignored in all of the paths above.
  - Reproducible runtime evidence:
    - I ran `./build/deterministic/tools/solver_cavity` on a temporary copy of `benchmarks/lid_driven_cavity_smoke.cfg` with `poisson_max_iterations = 1`.
    - The process exited successfully and advanced all 12 steps.
    - It reported:
      - `pressure_iterations: 1`
      - `pressure_relative_residual: 0.317217`
      - `divergence_l2: 0.00100295`
    - The configured pressure tolerance in that repro remained `1e-10`, so the driver demonstrably continued after a badly non-converged projection.

### 2. Medium: Taylor-Green temporal validation can silently underreport mismatched fields

- **Affected file / functions:**
  - `validation/run_validation_suite.py` - `read_vtk_velocity(...)`
  - `validation/run_validation_suite.py` - `relative_l2_difference(...)`
  - `validation/run_validation_suite.py` - `write_taylor_temporal_reports(...)`
- **Why this is wrong:**
  - `relative_l2_difference(left, right)` computes the numerator with `zip(left, right)` but never checks that the two vectors have the same length.
  - `read_vtk_velocity(...)` parses whatever velocity rows it finds in the VTK file and does not validate the expected vector count against the configured grid.
  - If one VTK is truncated or partially parsed, the temporal self-difference can be computed only on the shared prefix and still return a small or even zero value.
- **Impact:**
  - The temporal self-convergence gate can falsely pass even when one of the compared final fields is incomplete.
  - That undermines the correctness of the validation harness's reported temporal order.
- **Concrete evidence:**
  - Local repro:
    - `relative_l2_difference([1.0, 1000.0], [1.0])` currently returns `0.0`.
  - That result is mathematically wrong for mismatched fields; the extra `1000.0` entry in the left-hand vector is silently ignored.
  - `write_taylor_temporal_reports(...)` uses this helper directly for the Taylor-Green self-convergence table, so the bad comparison feeds a top-level validation gate.

## Notes

- `ctest --preset deterministic` passed during this review.
- The passing test suite does not invalidate the findings above; it shows these failure modes are not currently covered by the existing deterministic runtime tests.

## No Remediation Applied

This document records findings only. No code, configuration, thresholds, workflows, or generated artifacts were changed as part of this review beyond adding this Markdown file.
