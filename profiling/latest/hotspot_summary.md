# Default-Policy Hotspots

| Symbol | Samples | Share | Category |
| --- | --- | --- | --- |
| solver::linsolve::(anonymous namespace)::axpy_active(solver::StructuredField&, solver::StructuredField const&, double) | 12 | 0.078 | other |
| solver::(anonymous namespace)::solve_tridiagonal(std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>>&) | 9 | 0.058 | predictor_adi |
| solver::(anonymous namespace)::apply_face_field_boundary(solver::FaceField&, solver::BoundaryFace, solver::BoundaryCondition const&) | 8 | 0.052 | other |
| solver::linsolve::(anonymous namespace)::neighbor_value(solver::PressureField const&, solver::PressureBoundarySet const&, int, int, int, solver::Axis, bool) | 8 | 0.052 | pressure_solve |
| solver::(anonymous namespace)::reconstruct_face_y(solver::StructuredField const&, int, int, int, double, solver::AdvectionOptions const&) | 8 | 0.052 | advection |
| solver::(anonymous namespace)::reconstruct_face_x(solver::StructuredField const&, int, int, int, double, solver::AdvectionOptions const&) | 6 | 0.039 | advection |
| solver::linsolve::(anonymous namespace)::apply_poisson_operator(solver::PressureField const&, solver::PressureBoundarySet const&, solver::PressureField&) | 6 | 0.039 | pressure_solve |
| solver::linsolve::solve_pressure_poisson(solver::ScalarField const&, solver::PressureBoundarySet const&, solver::ProjectionOptions const&, solver::PressureField&) | 6 | 0.039 | pressure_solve |
| solver::linsolve::(anonymous namespace)::smooth_jacobi(solver::linsolve::(anonymous namespace)::HierarchyLevel&, solver::PressureBoundarySet const&, solver::linsolve::MultigridPolicy const&, int, bool) | 6 | 0.039 | pressure_solve |
| solver::(anonymous namespace)::compute_w_advection_2d(solver::VelocityField const&, solver::AdvectionOptions const&, solver::FaceField&) | 6 | 0.039 | advection |
| solver::linsolve::(anonymous namespace)::fill_storage(solver::StructuredField&, double) | 6 | 0.039 | other |
| solver::AlignedBuffer<double, 64ul>::reset() | 5 | 0.032 | other |
