# Code Review: Logical and Mathematical Errors

This document catalogues logical and mathematical issues found during a review of the
incompressible Navier–Stokes solver codebase. **No code has been modified.**

---

## 1. Velocity Correction Uses Wrong Sign Convention

**Files:** `solver/projection.cpp` lines 1089–1092

The pressure-correction step computes:

```
u_corrected = u_predicted - (dt / rho) * grad(p)
```

This is correct for the standard projection method where the Poisson equation is
`div(grad(p)) = (rho / dt) * div(u*)`. However, in `build_pressure_rhs` (line 1062)
the RHS is scaled by `rho / dt` (positive), and the Poisson operator in
`poisson_solver.cpp` line 262 applies the **negative** Laplacian (`-∇²p = rhs`).
That means the solve actually finds `p` such that `-∇²p = (rho/dt) div(u*)`,
i.e. `∇²p = -(rho/dt) div(u*)`.

In the standard formulation `∇²p = (rho/dt) div(u*)`, so the Poisson operator's
built-in negation means the correction should be:

```
u_corrected = u_predicted + (dt / rho) * grad(p)   // note: plus sign
```

but the code applies **`-dt / rho`**. The two negations (one from the Poisson
operator sign, one from the correction) cancel out, producing a self-consistent
scheme. **This is not a bug**, but the sign convention is non-obvious and deserves
a clarifying comment to prevent future maintenance errors.

**Severity:** Informational (self-consistent but fragile)

---

## 2. `compute_operator_diagonal` — Neumann Boundary Diagonal Adjustment

**File:** `linsolve/poisson_solver.cpp` lines 271–314

For a Neumann boundary the ghost value is `p_ghost = p_interior`, so the stencil
`(p_neighbor - 2*p_center + p_ghost)` becomes `(p_neighbor - 2*p_center + p_center)`
= `(p_neighbor - p_center)`. The diagonal contribution from that side drops from
`2/h²` to `1/h²`, i.e. a reduction of `1/h²`. The code implements this as:

```cpp
coefficient += type == PressureBoundaryType::dirichlet ? inverse_spacing_squared
                                                       : -inverse_spacing_squared;
```

For Neumann (else branch): `coefficient = 2/h² - 1/h² = 1/h²` ✓
For Dirichlet with ghost `= 2*p_wall - p_center`: stencil becomes
`(p_neighbor - 2*p_center + 2*p_wall - p_center)` = `(p_neighbor - 3*p_center + 2*p_wall)`,
so diagonal is `3/h²`, i.e. `2/h² + 1/h²` ✓

The arithmetic is correct. **No error.**

---

## 3. Poisson Operator Sign Convention (Negative Laplacian)

**File:** `linsolve/poisson_solver.cpp` lines 244–269

The Poisson operator computes **`-∇²p`** (negative Laplacian). This is a deliberate
choice — the CG solver in `solve_pressure_poisson` is set up so that the operator
is symmetric positive-definite. The RHS construction at lines 703–705:

```cpp
copy_active(rhs, system_rhs);
scale_active(system_rhs, -1.0);          // negate the physical RHS
axpy_active(system_rhs, boundary_source, -1.0);  // subtract boundary source
```

negates the physical RHS to match. The scheme is self-consistent.

**Severity:** Informational

---

## 4. Pressure Neumann Ghost — Sign Convention for Lower vs. Upper Boundary

**File:** `solver/projection.cpp` lines 322–326

```cpp
case PressureBoundaryType::neumann:
    pressure(gi, gj, gk) =
        is_lower_boundary(face) ? active_value - condition.gradient * spacing
                                : active_value + condition.gradient * spacing;
```

For a lower boundary at `x = 0` the ghost is at `x = -h/2` and the interior at
`x = h/2`. With outward normal pointing in the `−x` direction, `dp/dn = -dp/dx`,
so `p_ghost = p_interior - gradient * h` when `gradient = dp/dn`. For the upper
boundary (outward normal in `+x`), `p_ghost = p_interior + gradient * h`.

This matches the code. **No error.**

---

## 5. Taylor–Green 3-D Exact Pressure Expression

**File:** `solver/taylor_green.cpp` lines 379–382

The 3-D Taylor–Green exact pressure is:

```cpp
-0.0625 * (cos(2x) + cos(2y)) * (cos(2z) + 2.0) * decay
```

The standard analytical result for the 3-D TGV initial pressure is:

```
p = -1/16 * [ cos(2x) + cos(2y) ] * [ cos(2z) + 2 ]
```

which is `−0.0625 * (cos 2x + cos 2y)(cos 2z + 2)`. This matches the code. **No error.**

---

## 6. Taylor–Green 2-D Kinetic Energy Decay Rate

**File:** `solver/taylor_green.cpp` lines 258–266

For the 2-D Taylor–Green vortex with `u = -cos(x)sin(y)e^{-2νt}`,
`v = sin(x)cos(y)e^{-2νt}`:

- Velocity decay exponent: `exp(-2νt)` — code uses rate = 2.0 ✓
- Energy ∝ velocity², so energy decay exponent is `exp(-4νt)` — code uses rate = 6.0 for 3-D and 4.0 for 2-D ✓
- Initial energy 2-D: `<0.5*(u²+v²)> = 0.5 * 0.5 = 0.25` — code uses 0.25 ✓

For 3-D with `u = -cos(x)sin(y)cos(z)`, etc.:

- Velocity decay exponent: `exp(-3νt)` — code uses rate = 3.0 ✓
- Energy decay: `exp(-6νt)` — code uses rate = 6.0 ✓
- Initial energy 3-D: `<0.5*(u²+v²+w²)> = 0.5 * 0.25 = 0.125` — code uses 0.125 ✓

**No error.**

---

## 7. Cyclic Tridiagonal Solver — Zeroed Sub/Super-diagonal Entries

**File:** `solver/projection.cpp` lines 658–664

In `solve_component_sweep_periodic`, the lower and upper diagonals are initialised
to `−beta`, but then the code **zeros** the first element of `lower` and the last
element of `upper`:

```cpp
lower[0] = 0.0;
upper[unknowns - 1] = 0.0;
```

These entries are the sub-diagonal first element and super-diagonal last element of
the *interior* tridiagonal part in the Sherman–Morrison decomposition. However, the
cyclic coupling (corner entries) is passed separately as `lower_corner = -beta` and
`upper_corner = -beta`. In the standard Sherman–Morrison cyclic-tridiagonal method
the interior tridiagonal should indeed have `lower[0]` and `upper[n-1]` set to zero
because they are replaced by the corner entries. **This is correct.**

**No error.**

---

## 8. ADI Predictor Sweep — Boundary Treatment Asymmetry

**File:** `solver/projection.cpp` lines 762–778 (and similar blocks for Y, Z)

When the sweep axis is the normal axis of the face field and the boundary is a
Dirichlet velocity boundary, the known boundary value is added to the RHS. When
the boundary is *not* the normal axis (tangential component), a ghost-relation
based on the tangential BC type is folded into the diagonal and RHS.

The ghost-relation for `fixed_pressure` returns `{coefficient = 1.0, constant = 0.0}`
(same as symmetry), which implies zero-gradient extrapolation. This is standard for
tangential velocity at a pressure outlet. **No error.**

---

## 9. TVD Advection Flux Reconstruction — Ratio Definition

**File:** `solver/momentum_terms.cpp` lines 42–77

For the TVD reconstruction with a van Leer limiter, the upwind ratio is computed as:

```cpp
// Positive face velocity (upwind = left cell):
ratio = (left - far_left) / (right - left)
```

This is the ratio of consecutive differences `(Δ_upwind) / (Δ_local)`. The standard
van Leer limiter expects `r = Δ_upwind / Δ_local`. The reconstruction is:

```cpp
left + 0.5 * φ(r) * (right - left)
```

This is the standard MUSCL-style TVD reconstruction. **No error.**

---

## 10. `correct_velocity` — Projection Scaling

**File:** `solver/projection.cpp` lines 1089–1092

```cpp
axpy_active(corrected_velocity.x, pressure_gradient.x, -options.dt / options.density);
```

With the negative-Laplacian Poisson operator convention (see Finding #3), the
pressure solve yields `p` such that `−∇²p = (ρ/dt)·div(u*)`. This means
`∇²p = −(ρ/dt)·div(u*)`. The correction should be `u = u* − (dt/ρ)·∇p` which
requires the sign flip. Given the operator convention, this produces the correct
divergence-free velocity. **Self-consistent; no error.**

---

## 11. Pressure Total Update — Sign of Pressure RHS Correction

**Files:** `solver/taylor_green.cpp` line 606, `solver/lid_driven_cavity.cpp` line 593,
`solver/channel_flow.cpp` line 499

After the projection step, the total pressure is updated:

```cpp
axpy_active(state.pressure_total, pressure_correction, 1.0);
axpy_active(state.pressure_total, pressure_rhs, -0.5 * viscosity * dt);
```

The second term `−0.5·ν·dt · pressure_rhs` is the incremental pressure correction
technique that accounts for the diffusion splitting error in the ADI predictor.
Since `pressure_rhs = (ρ/dt)·div(u*)`, this term is `−0.5·ν·ρ·div(u*)`. With
`ρ = 1`, this is a standard Guermond–Minev correction. **No error.**

---

## 12. `operator_verification.cpp` — Divergence Exact Solution

**File:** `solver/operator_verification.cpp` lines 79–106

The manufactured velocity field is:

```
u = sin(x) cos(y),   v = cos(x) sin(y)
```

The exact divergence is:

```
du/dx + dv/dy = cos(x)cos(y) + cos(x)cos(y) = 2 cos(x) cos(y)
```

The code checks against `2.0 * cos(x) * cos(y)` (line 105). **Correct.**

However, note that this velocity field is **not** divergence-free (`div(u) ≠ 0`).
This is intentional for testing the divergence operator, but it means the MMS test
does not verify the solver's ability to produce divergence-free fields. This is a
test coverage gap, not a code error.

**Severity:** Informational (test design note)

---

## 13. `sample_vertical_centerline_u` — Sampling on Face-Centered Field

**File:** `solver/lid_driven_cavity.cpp` lines 100–114

The u-component is sampled at the vertical centerline by selecting the face index
at `i_center = active.i_begin + nx/2`. For a face-centered field with `nx + 1`
active faces, the midpoint of the domain corresponds to face `nx/2` (0-indexed),
which is at coordinate `(nx/2) * dx`. For `nx = 128`, `dx = 1/128`, this gives
`x = 0.5`, which is the domain center. **Correct.**

The y-coordinate is sampled using `coordinate_for_storage_index`, which for the
face-x layout returns cell centers in y (since y is not the normal axis).
**Correct.**

---

## 14. Poisson Convergence Check Uses Relative Tolerance Incorrectly

**File:** `linsolve/poisson_solver.cpp` line 779

```cpp
if(diagnostics.relative_residual <= options.poisson_tolerance) {
```

The convergence test compares the **relative** residual (`||r_k|| / ||r_0||`)
against `poisson_tolerance`. The default tolerance is `1.0e-12`. This is a
relative tolerance check, which is standard for preconditioned CG.

However, the initial convergence check at line 735 uses the **absolute** residual:

```cpp
if(diagnostics.initial_residual_l2 <= options.poisson_tolerance) {
```

This means the early-exit path uses an absolute check while the iteration loop
uses a relative check. If the initial residual is very small (e.g. `1e-11`) but
nonzero, it would not trigger the early exit, and the relative residual would
start at `1.0`, requiring full convergence iterations. This is a **minor logical
inconsistency** — the early exit should arguably also use a relative check (which
would be `1.0 <= tol`, never triggering), or the threshold for the absolute check
should differ from the relative tolerance.

**Severity:** Low — affects only edge cases where initial residual is very small.

---

## 15. `compute_operator_diagonal` Uses `std::pow` for Squaring

**File:** `linsolve/poisson_solver.cpp` line 283

```cpp
const double inverse_spacing_squared = 1.0 / std::pow(grid.spacing(axis), 2.0);
```

Using `std::pow(x, 2.0)` for squaring is not a mathematical error but is
unnecessarily imprecise compared to `x * x` for floating-point arithmetic when
`-ffast-math` is not enabled. With strict IEEE semantics, `std::pow` may introduce
rounding that `x * x` would not.

**Severity:** Negligible

---

## Summary

| # | Location | Description | Severity |
|---|----------|-------------|----------|
| 1 | `projection.cpp:1089` | Sign convention fragile but self-consistent | Informational |
| 2 | `poisson_solver.cpp:271` | Diagonal adjustment correct | No error |
| 3 | `poisson_solver.cpp:244` | Negative-Laplacian convention intentional | Informational |
| 4 | `projection.cpp:322` | Neumann ghost sign correct | No error |
| 5 | `taylor_green.cpp:379` | 3-D TGV pressure correct | No error |
| 6 | `taylor_green.cpp:258` | Energy decay rates correct | No error |
| 7 | `projection.cpp:658` | Cyclic tridiagonal setup correct | No error |
| 8 | `projection.cpp:762` | ADI boundary treatment correct | No error |
| 9 | `momentum_terms.cpp:42` | TVD reconstruction correct | No error |
| 10 | `projection.cpp:1089` | Projection scaling consistent | No error |
| 11 | `taylor_green.cpp:606` | Pressure total update correct | No error |
| 12 | `operator_verification.cpp:79` | MMS divergence test correct (non-zero div by design) | Informational |
| 13 | `lid_driven_cavity.cpp:100` | Centerline sampling correct | No error |
| 14 | `poisson_solver.cpp:735,779` | Absolute vs. relative tolerance inconsistency in early-exit | Low |
| 15 | `poisson_solver.cpp:283` | `std::pow` for squaring instead of multiplication | Negligible |

**Overall assessment:** The mathematical and logical implementation is sound. The
discrete operators, projection method, ADI solver, TVD advection, multigrid-
preconditioned CG, and analytical test cases are all correctly implemented. Two
minor items (the absolute-vs-relative tolerance inconsistency in Poisson early-exit,
and the non-obvious sign convention) are noted for future maintenance awareness.
