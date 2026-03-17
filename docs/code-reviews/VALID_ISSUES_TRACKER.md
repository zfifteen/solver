# Valid Issues Tracker

**Date:** 2026-03-17  
**Purpose:** Standalone actionable-issues record derived from the review documents under `docs/code-reviews`. This tracker is intended to preserve the triage outcome and enough source-detail context that the individual review documents can be deleted without losing decision history.

## Triage Policy

- **Included:** confirmed actionable defects that can change solver correctness, convergence behavior, or validation outcomes.
- **Excluded:** `no error` findings, informational notes, design tradeoffs, optimization ideas, verification suggestions, naming/documentation improvements, and low-severity maintenance inconsistencies.
- **Deduplication rule:** if multiple review documents report the same valid issue, the tracker keeps one canonical issue entry and records the duplicate source references under it.

## Source Review Verdicts

### `CODE_REVIEW_FINDINGS.md`

- **Document type:** mixed mathematical review, mostly “no error” validation with a few informational notes.
- **Accepted actionable defects:** 0
- **Rejected items:** 15
- **Overall verdict:** no entries from this document qualify for the actionable tracker.

#### Item-by-item disposition

| # | Source finding | Outcome | Reason |
| --- | --- | --- | --- |
| 1 | Velocity correction uses wrong sign convention | Rejected | Informational only; the document itself concludes the scheme is self-consistent, not incorrect. |
| 2 | `compute_operator_diagonal` Neumann boundary adjustment | Rejected | No bug; the document concludes the arithmetic is correct. |
| 3 | Poisson operator sign convention | Rejected | Informational only; deliberate negative-Laplacian SPD formulation. |
| 4 | Pressure Neumann ghost sign convention | Rejected | No bug; document concludes the implementation is correct. |
| 5 | Taylor-Green 3D exact pressure expression | Rejected | No bug. |
| 6 | Taylor-Green 2D/3D energy decay rates | Rejected | No bug. |
| 7 | Cyclic tridiagonal zeroed sub/super-diagonal entries | Rejected | No bug; correct Sherman-Morrison setup. |
| 8 | ADI predictor boundary treatment asymmetry | Rejected | No bug; document concludes the behavior is standard and correct. |
| 9 | TVD advection flux ratio definition | Rejected | No bug; standard MUSCL/TVD ratio. |
| 10 | `correct_velocity` projection scaling | Rejected | No bug; same sign-convention note as item 1. |
| 11 | Pressure-total update sign of pressure-RHS term | Rejected | No bug; document classifies it as standard incremental-pressure correction. |
| 12 | Divergence MMS exact solution is nonzero | Rejected | Informational only; intentional test design note, not solver defect. |
| 13 | Cavity centerline sampling on face-centered field | Rejected | No bug. |
| 14 | Poisson convergence check uses absolute early-exit then relative iterative tolerance | Rejected | Real but low-severity maintenance inconsistency; excluded by cutoff because it does not create a demonstrated correctness/convergence failure in supported runs. |
| 15 | `std::pow` used for squaring | Rejected | Negligible implementation style/perf note, not a defect. |

### `LOGICAL_MATHEMATICAL_CODE_REVIEW.md`

- **Document type:** findings-focused confirmed-issues review.
- **Accepted actionable defects:** 2
- **Rejected items:** 0
- **Overall verdict:** both recorded findings are valid and retained.

#### Accepted findings from this source

| Source finding | Canonical tracker issue |
| --- | --- |
| CPU driver loops commit states after non-converged pressure projections | Issue 1 |
| Taylor-Green temporal validation can silently underreport mismatched fields | Issue 2 |

### `LOGICAL_MATHEMATICAL_CODE_REVIEW_GROK.md`

- **Document type:** findings-focused confirmed-issues review with one extra backend concern.
- **Accepted actionable defects:** 2 duplicates
- **Rejected items:** 1
- **Overall verdict:** two findings duplicate the accepted issues above; one item is excluded.

#### Item-by-item disposition

| Source finding | Outcome | Reason |
| --- | --- | --- |
| Drivers commit state after non-converged pressure projections | Accepted as duplicate | Same defect as Issue 1. |
| Validation suite L2 difference computation is unsafe | Accepted as duplicate | Same defect as Issue 2. |
| Metal backend precision/tolerance inconsistency | Rejected | Backend tradeoff for the intentionally narrow supported Metal path; not a confirmed correctness defect in accepted usage. |

### `poisson-solver-review.md`

- **Document type:** architectural assessment and recommendation memo for the Poisson solver.
- **Accepted actionable defects:** 0
- **Rejected items:** review-wide recommendation set
- **Overall verdict:** no confirmed solver defects to track; the document is preserved in this tracker as recommendation context only.

#### Review-wide disposition

The document praises the solver as numerically sound and “production-ready,” then proposes follow-up work such as:

- verifying loop stride against storage layout
- profiling the coarse solve contribution
- adding an operator symmetry verification test
- tracking multigrid iteration count versus grid size
- documenting coarsest-grid policy
- experimenting with smoother tuning, SIMD hints, and alternate prolongation

These are verification, profiling, tuning, or documentation suggestions, not confirmed defects, so none are retained as actionable issues.

## Valid Issues

### Issue 1: Drivers commit state after non-converged pressure projections

- **Severity:** High
- **Source documents:** `LOGICAL_MATHEMATICAL_CODE_REVIEW.md`, `LOGICAL_MATHEMATICAL_CODE_REVIEW_GROK.md`
- **Affected files / behavior:**
  - `solver/channel_flow.cpp`: timestep loop continues after `project_velocity(...)`
  - `solver/lid_driven_cavity.cpp`: timestep loop continues after `project_velocity(...)`
  - `solver/taylor_green.cpp`: CPU timestep loop continues after `project_velocity(...)`
  - `solver/taylor_green.cpp`: final CPU cleanup in the Metal path also updates final state without checking convergence
- **Why this is a valid defect:**
  - `project_velocity(...)` returns `ProjectionDiagnostics`, including `pressure_solve.converged`, `pressure_solve.iterations`, and `pressure_solve.relative_residual`.
  - The callers consume the corrected velocity and update total pressure immediately, but use the returned pressure diagnostics only for reporting.
  - If the pressure solve fails to meet tolerance, the driver still commits the updated state and keeps advancing, which violates the intended projection contract and can leave the state materially under-projected.
- **Concrete evidence retained from review work:**
  - Code inspection confirmed the convergence flag is ignored in the driver paths listed above.
  - Reproduced behavior on cavity smoke run using a temporary config with `poisson_max_iterations = 1`:
    - run completed successfully
    - `pressure_iterations: 1`
    - `pressure_relative_residual: 0.317217`
    - `divergence_l2: 0.00100295`
    - configured pressure tolerance remained `1e-10`
  - This demonstrates that a badly non-converged projection can be committed as if the run succeeded normally.

### Issue 2: Taylor-Green temporal validation can silently underreport mismatched fields

- **Severity:** Medium
- **Source documents:** `LOGICAL_MATHEMATICAL_CODE_REVIEW.md`, `LOGICAL_MATHEMATICAL_CODE_REVIEW_GROK.md`
- **Affected files / behavior:**
  - `validation/run_validation_suite.py`
  - functions involved: `read_vtk_velocity(...)`, `relative_l2_difference(...)`, `write_taylor_temporal_reports(...)`
- **Why this is a valid defect:**
  - `relative_l2_difference(left, right)` computes the numerator with `zip(left, right)` and does not check that the vectors have the same length.
  - `read_vtk_velocity(...)` collects whatever velocity rows it finds and does not validate the parsed vector count against the expected field size.
  - If one VTK output is truncated or partially parsed, the comparison can silently operate on only the common prefix and still produce an artificially small or zero difference.
- **Concrete evidence retained from review work:**
  - Direct repro of the current helper behavior:
    - `relative_l2_difference([1.0, 1000.0], [1.0]) == 0.0`
  - That result is mathematically wrong for mismatched fields because the extra left-hand entry is ignored entirely.
  - `write_taylor_temporal_reports(...)` uses this helper directly to compute the Taylor-Green temporal self-convergence gate, so the defect can produce false-positive validation outcomes.

## Excluded Items

The following categories were reviewed and intentionally excluded from the actionable tracker:

- **Self-consistent sign-convention notes**
  - pressure/projection sign convention
  - negative-Laplacian Poisson convention
  - projection scaling note
- **Findings explicitly concluded as correct**
  - Neumann ghost handling
  - operator diagonal adjustments
  - Taylor-Green exact pressure/energy formulas
  - cyclic tridiagonal periodic setup
  - TVD limiter ratio definition
  - cavity centerline sampling
- **Test-design or maintenance notes**
  - divergence MMS uses a non-divergence-free manufactured field by design
  - Poisson absolute-versus-relative tolerance inconsistency on early exit
  - `std::pow` used for squaring
- **Backend tradeoffs / supported-path limitations**
  - Metal float/double precision mismatch concern
  - Metal tolerance floor concerns in the narrow periodic 3D Taylor-Green path
- **Recommendations rather than defects**
  - profiling actions
  - symmetry/verification additions
  - multigrid iteration tracking
  - coarsening policy documentation
  - smoother/prolongation/SIMD tuning ideas

## Deletion Readiness

This tracker now captures:

- which review documents were examined
- the disposition of each document
- the disposition of every named finding from the defect-oriented review docs
- the two canonical actionable issues that remain after deduplication
- the reason the rest of the material was excluded

With this information preserved here, the individual review documents can be deleted without losing the actionable defect set or the rationale for rejecting the other review claims.

## No Remediation Applied

This tracker is a triage artifact only. No solver code, tests, configs, workflows, or source review documents were modified while preparing it, beyond updating this tracker file itself.
