# Default-Policy Hotspots

| Symbol | Samples | Share | Category |
| --- | --- | --- | --- |
| solver::(anonymous namespace)::solve_tridiagonal(std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>>&) | 11 | 0.061 | predictor_adi |
| solver::(anonymous namespace)::apply_face_field_boundary(solver::FaceField&, solver::BoundaryFace, solver::BoundaryCondition const&) | 10 | 0.055 | other |
| solver::linsolve::(anonymous namespace)::neighbor_value(solver::PressureField const&, solver::PressureBoundarySet const&, int, int, int, solver::Axis, bool) | 8 | 0.044 | pressure_solve |
| solver::linsolve::(anonymous namespace)::axpy_active(solver::StructuredField&, solver::StructuredField const&, double) | 8 | 0.044 | other |
| solver::(anonymous namespace)::compute_u_advection(solver::VelocityField const&, solver::AdvectionOptions const&, solver::FaceField&) | 6 | 0.033 | other |
| solver::(anonymous namespace)::solve_component_sweep(solver::FaceField const&, solver::Axis, solver::BoundaryConditionSet const&, double, solver::FaceField&, int&) | 6 | 0.033 | predictor_adi |
| solver::linsolve::(anonymous namespace)::active_mean(solver::StructuredField const&) | 6 | 0.033 | other |
| solver::linsolve::(anonymous namespace)::smooth_jacobi(solver::linsolve::(anonymous namespace)::HierarchyLevel&, solver::PressureBoundarySet const&, solver::linsolve::MultigridPolicy const&, int, bool) | 6 | 0.033 | pressure_solve |
| solver::linsolve::solve_pressure_poisson(solver::ScalarField const&, solver::PressureBoundarySet const&, solver::ProjectionOptions const&, solver::PressureField&) | 6 | 0.033 | pressure_solve |
| solver::linsolve::(anonymous namespace)::apply_poisson_operator(solver::PressureField const&, solver::PressureBoundarySet const&, solver::PressureField&) | 5 | 0.028 | pressure_solve |
| solver::(anonymous namespace)::compute_w_advection(solver::VelocityField const&, solver::AdvectionOptions const&, solver::FaceField&) | 5 | 0.028 | other |
| solver::linsolve::(anonymous namespace)::prolongate_and_add(solver::PressureField const&, solver::PressureField&) | 5 | 0.028 | other |
