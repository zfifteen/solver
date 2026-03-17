# Default-Policy Hotspots

| Symbol | Samples | Share | Category |
| --- | --- | --- | --- |
| solver::linsolve::(anonymous namespace)::axpy_active(solver::StructuredField&, solver::StructuredField const&, double) | 12 | 0.065 | other |
| solver::(anonymous namespace)::apply_face_field_boundary(solver::FaceField&, solver::BoundaryFace, solver::BoundaryCondition const&) | 11 | 0.059 | other |
| solver::(anonymous namespace)::solve_tridiagonal(std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>> const&, std::__1::vector<double, std::__1::allocator<double>>&) | 11 | 0.059 | predictor_adi |
| solver::linsolve::(anonymous namespace)::smooth_jacobi(solver::linsolve::(anonymous namespace)::HierarchyLevel&, solver::PressureBoundarySet const&, solver::linsolve::MultigridPolicy const&, int, bool) | 9 | 0.048 | pressure_solve |
| solver::linsolve::(anonymous namespace)::neighbor_value(solver::PressureField const&, solver::PressureBoundarySet const&, int, int, int, solver::Axis, bool) | 8 | 0.043 | pressure_solve |
| solver::linsolve::(anonymous namespace)::apply_poisson_operator(solver::PressureField const&, solver::PressureBoundarySet const&, solver::PressureField&) | 6 | 0.032 | pressure_solve |
| solver::apply_velocity_boundary_conditions(solver::BoundaryConditionSet const&, solver::VelocityField&) | 6 | 0.032 | other |
| solver::AlignedBuffer<double, 64ul>::resize(unsigned long) | 5 | 0.027 | other |
| solver::(anonymous namespace)::axpy_active(solver::StructuredField&, solver::StructuredField const&, double) | 5 | 0.027 | other |
| solver::linsolve::solve_pressure_poisson(solver::ScalarField const&, solver::PressureBoundarySet const&, solver::ProjectionOptions const&, solver::PressureField&) | 5 | 0.027 | pressure_solve |
| solver::linsolve::(anonymous namespace)::active_mean(solver::StructuredField const&) | 5 | 0.027 | other |
| solver::(anonymous namespace)::compute_u_advection(solver::VelocityField const&, solver::AdvectionOptions const&, solver::FaceField&) | 5 | 0.027 | other |
