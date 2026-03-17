# Fix Plan: Taylor-Green Validation Must Reject Mismatched Field Lengths

**Related tracker issue:** Issue 2 in [VALID_ISSUES_TRACKER.md](/Users/velocityworks/IdeaProjects/solver/docs/code-reviews/VALID_ISSUES_TRACKER.md)

## Summary

Make the Taylor-Green temporal self-convergence path fail closed when a VTK field is truncated, partially parsed, or length-mismatched. The validation harness should raise a clear error instead of silently comparing only the shared prefix of two vectors.

## Implementation Changes

### 1. Make the vector-difference helper strict

- Update `relative_l2_difference(...)` in [run_validation_suite.py](/Users/velocityworks/IdeaProjects/solver/validation/run_validation_suite.py) to validate:
  - both inputs are non-empty
  - both inputs have equal length
- On violation, raise `ValueError` with a message that states the two lengths or identifies the empty-input condition.
- Keep the actual L2 formula unchanged for valid inputs.

### 2. Validate parsed VTK field size

- Change `read_vtk_velocity(...)` so it accepts the expected velocity field size for the case being read.
- In the Taylor-Green temporal study, derive the expected scalar-component count from the generated case dimensions:
  - expected velocity vectors = `nx * ny * nz`
  - expected flat list length = `3 * nx * ny * nz`
- After parsing, verify that the flattened velocity list exactly matches that expected length.
- If the count is short, long, or missing entirely, raise `ValueError` with the file path and both expected/actual counts.

### 3. Thread the stricter contract through the temporal validation path

- Update `write_taylor_temporal_reports(...)` to pass the expected count into `read_vtk_velocity(...)`.
- Allow malformed VTK data to abort the validation suite rather than writing a misleading temporal-convergence report.
- Keep all current report formats and thresholds unchanged for valid inputs.

## Interfaces and Behavior

- No solver code changes.
- No benchmark-threshold changes.
- Internal validation-script behavior changes:
  - previously, malformed/truncated velocity fields could produce false-small self-differences
  - after the fix, malformed/truncated velocity fields must raise and fail the validation run
- Generated report file formats remain unchanged on successful runs.

## Test Plan

- Add focused Python-level coverage for the validation helpers, either as a small `python3`-driven test script or a `ctest`-invoked Python test, covering:
  - equal-length vectors with known result
  - mismatched vector lengths raising `ValueError`
  - empty vectors raising `ValueError`
  - a truncated VTK fixture or synthetic sample raising on wrong parsed count
- Re-run the validation harness on a normal deterministic build and confirm:
  - the temporal report still generates successfully on valid data
  - the summary/gate behavior is unchanged for good inputs
- Acceptance criteria:
  - malformed or mismatched field data can no longer yield a numeric self-difference
  - valid Taylor-Green temporal validation still completes and preserves current thresholds

## Assumptions

- The correct behavior for malformed validation inputs is hard failure, not best-effort comparison.
- The expected field size for the temporal study should be derived from the generated config, not inferred from the file contents alone.
- This plan intentionally keeps the change narrow; it does not redesign the report generator or broaden validation-suit architecture.
