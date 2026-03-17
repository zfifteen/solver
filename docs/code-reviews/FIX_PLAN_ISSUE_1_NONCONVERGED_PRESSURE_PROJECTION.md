# Fix Plan: Non-Converged Pressure Projections Committed by Drivers

**Related tracker issue:** Issue 1 in [VALID_ISSUES_TRACKER.md](/Users/velocityworks/IdeaProjects/solver/docs/code-reviews/VALID_ISSUES_TRACKER.md)

## Summary

Prevent the simulation drivers from committing velocity and pressure state after `project_velocity(...)` reports a non-converged pressure solve. The fix should make projection acceptance explicit, fail closed on bad projections, and revalidate all supported deterministic benchmark paths before landing.

## Implementation Changes

### 1. Add a shared projection-acceptance guard

- Add a small shared helper in the solver layer that checks `ProjectionDiagnostics.pressure_solve.converged`.
- The helper should throw `std::runtime_error` when the pressure solve did not converge.
- The error message should include:
  - solver/driver name
  - step index
  - iteration count
  - reported relative residual
  - configured pressure tolerance
- Do not silently retry or auto-relax tolerance inside the helper.

### 2. Enforce the guard before state commit

- Call the helper immediately after each `project_velocity(...)` call and before any of the following state mutations:
  - updating corrected velocity into the persistent state
  - applying `pressure_correction` into total pressure
  - applying the `pressure_rhs` correction term into total pressure
- Apply this in:
  - [channel_flow.cpp](/Users/velocityworks/IdeaProjects/solver/solver/channel_flow.cpp)
  - [lid_driven_cavity.cpp](/Users/velocityworks/IdeaProjects/solver/solver/lid_driven_cavity.cpp)
  - [taylor_green.cpp](/Users/velocityworks/IdeaProjects/solver/solver/taylor_green.cpp) CPU timestep loop
  - the final CPU cleanup projection inside the Metal path in [taylor_green.cpp](/Users/velocityworks/IdeaProjects/solver/solver/taylor_green.cpp)
- Keep the existing non-finite-state checks; this change is additive and specifically closes the “non-converged but still committed” hole.

### 3. Preserve supported benchmark behavior explicitly

- After the fail-closed behavior is added, rerun the deterministic benchmark surface.
- If any supported benchmark now fails only because its configured pressure iteration budget is too small, increase that benchmark/config iteration budget rather than weakening the new guard.
- Treat benchmark/config updates as acceptable only if:
  - the solver then converges cleanly under the requested tolerance
  - the benchmark still satisfies its existing validation thresholds
- Do not add a fallback mode that commits non-converged states.

## Interfaces and Behavior

- No public API additions are required.
- Runtime behavior changes:
  - previously, some runs could complete after a failed pressure projection
  - after the fix, those runs must fail with a clear error instead of committing the bad state
- Existing result structs and printed metrics can stay unchanged.

## Test Plan

- Add runtime tests in [test_runtime.cpp](/Users/velocityworks/IdeaProjects/solver/tests/test_runtime.cpp) that intentionally starve the pressure solve (`poisson_max_iterations = 1`) and assert rejection for:
  - channel flow
  - lid-driven cavity
  - CPU Taylor-Green
- Add one Metal-path regression that exercises the final CPU cleanup rejection path and asserts the Metal run fails rather than publishing a cleaned-up final state from a non-converged cleanup projection.
- Re-run:
  - `ctest --preset deterministic`
  - the manual deterministic validation harness
  - the supported benchmark configs used by the validation harness
- Acceptance criteria:
  - starved-pressure tests fail closed with clear error messages
  - supported reference benchmarks still pass after any required iteration-budget adjustments
  - no driver path commits a state after `pressure_solve.converged == false`

## Assumptions

- The correct long-term contract is “do not commit non-converged projections.”
- If a currently supported benchmark misses tolerance, the preferred remedy is to raise its iteration budget or otherwise make it truly converge, not to retain permissive commit behavior.
- This plan does not introduce retries, warnings-only behavior, or new convergence-policy knobs.
