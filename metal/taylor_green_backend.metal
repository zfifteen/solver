#include <metal_stdlib>

using namespace metal;

#if SOLVER_METAL_USE_FLOAT
using solver_scalar = float;
#else
using solver_scalar = double;
#endif

constant solver_scalar kSmallDenominator = static_cast<solver_scalar>(1.0e-14);
constant uint kReductionWidth = 256;

struct SolverParams {
  uint nx;
  uint ny;
  uint nz;
  uint use_ab2;
  uint advection_scheme;
  solver_scalar dx;
  solver_scalar dy;
  solver_scalar dz;
  solver_scalar viscosity;
  solver_scalar dt;
  solver_scalar alpha;
};

inline uint wrap_index(const int index, const uint extent) {
  const int wrapped = index % static_cast<int>(extent);
  return static_cast<uint>(wrapped < 0 ? wrapped + static_cast<int>(extent) : wrapped);
}

inline uint linear_index(const uint i, const uint j, const uint k, constant SolverParams& params) {
  return i + params.nx * (j + params.ny * k);
}

inline void unpack_index(const uint gid,
                         constant SolverParams& params,
                         thread uint& i,
                         thread uint& j,
                         thread uint& k) {
  const uint plane = params.nx * params.ny;
  k = gid / plane;
  const uint rem = gid - k * plane;
  j = rem / params.nx;
  i = rem - j * params.nx;
}

inline solver_scalar load_field(const device solver_scalar* field,
                                const int i,
                                const int j,
                                const int k,
                                constant SolverParams& params) {
  return field[linear_index(wrap_index(i, params.nx),
                            wrap_index(j, params.ny),
                            wrap_index(k, params.nz),
                            params)];
}

inline solver_scalar limiter_van_leer(const solver_scalar ratio) {
  return (ratio + fabs(ratio)) / (static_cast<solver_scalar>(1.0) + fabs(ratio));
}

inline solver_scalar safe_ratio(const solver_scalar numerator, const solver_scalar denominator) {
  if(fabs(denominator) < kSmallDenominator) {
    return static_cast<solver_scalar>(0.0);
  }
  return numerator / denominator;
}

inline solver_scalar reconstruct_face_x(const device solver_scalar* field,
                                        const int i_right,
                                        const uint j,
                                        const uint k,
                                        const solver_scalar face_velocity,
                                        constant SolverParams& params) {
  const solver_scalar left = load_field(field, i_right - 1, static_cast<int>(j), static_cast<int>(k), params);
  const solver_scalar right = load_field(field, i_right, static_cast<int>(j), static_cast<int>(k), params);

  if(params.advection_scheme == 2u) {
    return static_cast<solver_scalar>(0.5) * (left + right);
  }

  if(params.advection_scheme == 1u || fabs(face_velocity) < kSmallDenominator) {
    return face_velocity >= static_cast<solver_scalar>(0.0) ? left : right;
  }

  if(face_velocity >= static_cast<solver_scalar>(0.0)) {
    const solver_scalar delta = right - left;
    const solver_scalar ratio = safe_ratio(
        left - load_field(field, i_right - 2, static_cast<int>(j), static_cast<int>(k), params),
        delta);
    return left + static_cast<solver_scalar>(0.5) * limiter_van_leer(ratio) * delta;
  }

  const solver_scalar delta = right - left;
  const solver_scalar ratio = safe_ratio(
      load_field(field, i_right + 1, static_cast<int>(j), static_cast<int>(k), params) - right,
      delta);
  return right - static_cast<solver_scalar>(0.5) * limiter_van_leer(ratio) * delta;
}

inline solver_scalar reconstruct_face_y(const device solver_scalar* field,
                                        const uint i,
                                        const int j_right,
                                        const uint k,
                                        const solver_scalar face_velocity,
                                        constant SolverParams& params) {
  const solver_scalar lower = load_field(field, static_cast<int>(i), j_right - 1, static_cast<int>(k), params);
  const solver_scalar upper = load_field(field, static_cast<int>(i), j_right, static_cast<int>(k), params);

  if(params.advection_scheme == 2u) {
    return static_cast<solver_scalar>(0.5) * (lower + upper);
  }

  if(params.advection_scheme == 1u || fabs(face_velocity) < kSmallDenominator) {
    return face_velocity >= static_cast<solver_scalar>(0.0) ? lower : upper;
  }

  if(face_velocity >= static_cast<solver_scalar>(0.0)) {
    const solver_scalar delta = upper - lower;
    const solver_scalar ratio = safe_ratio(
        lower - load_field(field, static_cast<int>(i), j_right - 2, static_cast<int>(k), params),
        delta);
    return lower + static_cast<solver_scalar>(0.5) * limiter_van_leer(ratio) * delta;
  }

  const solver_scalar delta = upper - lower;
  const solver_scalar ratio = safe_ratio(
      load_field(field, static_cast<int>(i), j_right + 1, static_cast<int>(k), params) - upper,
      delta);
  return upper - static_cast<solver_scalar>(0.5) * limiter_van_leer(ratio) * delta;
}

inline solver_scalar reconstruct_face_z(const device solver_scalar* field,
                                        const uint i,
                                        const uint j,
                                        const int k_right,
                                        const solver_scalar face_velocity,
                                        constant SolverParams& params) {
  const solver_scalar back = load_field(field, static_cast<int>(i), static_cast<int>(j), k_right - 1, params);
  const solver_scalar front = load_field(field, static_cast<int>(i), static_cast<int>(j), k_right, params);

  if(params.advection_scheme == 2u) {
    return static_cast<solver_scalar>(0.5) * (back + front);
  }

  if(params.advection_scheme == 1u || fabs(face_velocity) < kSmallDenominator) {
    return face_velocity >= static_cast<solver_scalar>(0.0) ? back : front;
  }

  if(face_velocity >= static_cast<solver_scalar>(0.0)) {
    const solver_scalar delta = front - back;
    const solver_scalar ratio = safe_ratio(
        back - load_field(field, static_cast<int>(i), static_cast<int>(j), k_right - 2, params),
        delta);
    return back + static_cast<solver_scalar>(0.5) * limiter_van_leer(ratio) * delta;
  }

  const solver_scalar delta = front - back;
  const solver_scalar ratio = safe_ratio(
      load_field(field, static_cast<int>(i), static_cast<int>(j), k_right + 1, params) - front,
      delta);
  return front - static_cast<solver_scalar>(0.5) * limiter_van_leer(ratio) * delta;
}

inline solver_scalar second_derivative_x(const device solver_scalar* field,
                                         const int i,
                                         const int j,
                                         const int k,
                                         constant SolverParams& params) {
  const solver_scalar inverse_spacing_squared = static_cast<solver_scalar>(1.0) / (params.dx * params.dx);
  return (load_field(field, i + 1, j, k, params) -
          static_cast<solver_scalar>(2.0) * load_field(field, i, j, k, params) +
          load_field(field, i - 1, j, k, params)) *
         inverse_spacing_squared;
}

inline solver_scalar second_derivative_y(const device solver_scalar* field,
                                         const int i,
                                         const int j,
                                         const int k,
                                         constant SolverParams& params) {
  const solver_scalar inverse_spacing_squared = static_cast<solver_scalar>(1.0) / (params.dy * params.dy);
  return (load_field(field, i, j + 1, k, params) -
          static_cast<solver_scalar>(2.0) * load_field(field, i, j, k, params) +
          load_field(field, i, j - 1, k, params)) *
         inverse_spacing_squared;
}

inline solver_scalar second_derivative_z(const device solver_scalar* field,
                                         const int i,
                                         const int j,
                                         const int k,
                                         constant SolverParams& params) {
  const solver_scalar inverse_spacing_squared = static_cast<solver_scalar>(1.0) / (params.dz * params.dz);
  return (load_field(field, i, j, k + 1, params) -
          static_cast<solver_scalar>(2.0) * load_field(field, i, j, k, params) +
          load_field(field, i, j, k - 1, params)) *
         inverse_spacing_squared;
}

inline solver_scalar second_derivative_axis(const device solver_scalar* field,
                                            const int i,
                                            const int j,
                                            const int k,
                                            const uint axis,
                                            constant SolverParams& params) {
  switch(axis) {
    case 0u:
      return second_derivative_x(field, i, j, k, params);
    case 1u:
      return second_derivative_y(field, i, j, k, params);
    default:
      return second_derivative_z(field, i, j, k, params);
  }
}

inline solver_scalar mixed_second_derivative(const device solver_scalar* field,
                                             const int i,
                                             const int j,
                                             const int k,
                                             const uint outer_axis,
                                             const uint inner_axis,
                                             constant SolverParams& params) {
  const solver_scalar spacing = outer_axis == 0u ? params.dx : outer_axis == 1u ? params.dy : params.dz;
  const solver_scalar inverse_spacing_squared = static_cast<solver_scalar>(1.0) / (spacing * spacing);
  const int di = outer_axis == 0u ? 1 : 0;
  const int dj = outer_axis == 1u ? 1 : 0;
  const int dk = outer_axis == 2u ? 1 : 0;

  return (second_derivative_axis(field, i + di, j + dj, k + dk, inner_axis, params) -
          static_cast<solver_scalar>(2.0) * second_derivative_axis(field, i, j, k, inner_axis, params) +
          second_derivative_axis(field, i - di, j - dj, k - dk, inner_axis, params)) *
         inverse_spacing_squared;
}

inline solver_scalar triple_second_derivative(const device solver_scalar* field,
                                              const int i,
                                              const int j,
                                              const int k,
                                              const uint axis_a,
                                              const uint axis_b,
                                              const uint axis_c,
                                              constant SolverParams& params) {
  const solver_scalar spacing = axis_a == 0u ? params.dx : axis_a == 1u ? params.dy : params.dz;
  const solver_scalar inverse_spacing_squared = static_cast<solver_scalar>(1.0) / (spacing * spacing);
  const int di = axis_a == 0u ? 1 : 0;
  const int dj = axis_a == 1u ? 1 : 0;
  const int dk = axis_a == 2u ? 1 : 0;

  return (mixed_second_derivative(field, i + di, j + dj, k + dk, axis_b, axis_c, params) -
          static_cast<solver_scalar>(2.0) * mixed_second_derivative(field, i, j, k, axis_b, axis_c, params) +
          mixed_second_derivative(field, i - di, j - dj, k - dk, axis_b, axis_c, params)) *
         inverse_spacing_squared;
}

inline solver_scalar factorized_operator_value(const device solver_scalar* field,
                                               const int i,
                                               const int j,
                                               const int k,
                                               constant SolverParams& params) {
  const solver_scalar center = load_field(field, i, j, k, params);
  const solver_scalar laplacian = second_derivative_x(field, i, j, k, params) +
                                  second_derivative_y(field, i, j, k, params) +
                                  second_derivative_z(field, i, j, k, params);
  const solver_scalar mixed = mixed_second_derivative(field, i, j, k, 0u, 1u, params) +
                              mixed_second_derivative(field, i, j, k, 0u, 2u, params) +
                              mixed_second_derivative(field, i, j, k, 1u, 2u, params);
  const solver_scalar triple = triple_second_derivative(field, i, j, k, 0u, 1u, 2u, params);
  return center - params.alpha * laplacian +
         params.alpha * params.alpha * mixed -
         params.alpha * params.alpha * params.alpha * triple;
}

inline solver_scalar pressure_gradient_u(const device solver_scalar* pressure,
                                         const int i,
                                         const int j,
                                         const int k,
                                         constant SolverParams& params) {
  return (load_field(pressure, i, j, k, params) -
          load_field(pressure, i - 1, j, k, params)) / params.dx;
}

inline solver_scalar pressure_gradient_v(const device solver_scalar* pressure,
                                         const int i,
                                         const int j,
                                         const int k,
                                         constant SolverParams& params) {
  return (load_field(pressure, i, j, k, params) -
          load_field(pressure, i, j - 1, k, params)) / params.dy;
}

inline solver_scalar pressure_gradient_w(const device solver_scalar* pressure,
                                         const int i,
                                         const int j,
                                         const int k,
                                         constant SolverParams& params) {
  return (load_field(pressure, i, j, k, params) -
          load_field(pressure, i, j, k - 1, params)) / params.dz;
}

kernel void init_taylor_green_3d(device solver_scalar* u [[buffer(0)]],
                                 device solver_scalar* v [[buffer(1)]],
                                 device solver_scalar* w [[buffer(2)]],
                                 device solver_scalar* pressure [[buffer(3)]],
                                 constant SolverParams& params [[buffer(4)]],
                                 uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }

  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);

  const solver_scalar x_face = static_cast<solver_scalar>(i) * params.dx;
  const solver_scalar y_face = static_cast<solver_scalar>(j) * params.dy;
  const solver_scalar x_cell = (static_cast<solver_scalar>(i) + static_cast<solver_scalar>(0.5)) * params.dx;
  const solver_scalar y_cell = (static_cast<solver_scalar>(j) + static_cast<solver_scalar>(0.5)) * params.dy;
  const solver_scalar z_cell = (static_cast<solver_scalar>(k) + static_cast<solver_scalar>(0.5)) * params.dz;

  u[gid] = -cos(x_face) * sin(y_cell) * cos(z_cell);
  v[gid] = sin(x_cell) * cos(y_face) * cos(z_cell);
  w[gid] = static_cast<solver_scalar>(0.0);
  pressure[gid] =
      static_cast<solver_scalar>(-0.0625) *
      (cos(static_cast<solver_scalar>(2.0) * x_cell) + cos(static_cast<solver_scalar>(2.0) * y_cell)) *
      (cos(static_cast<solver_scalar>(2.0) * z_cell) + static_cast<solver_scalar>(2.0));
}

kernel void assemble_rhs_u(device solver_scalar* advection_out [[buffer(0)]],
                           device solver_scalar* rhs_out [[buffer(1)]],
                           const device solver_scalar* current_u [[buffer(2)]],
                           const device solver_scalar* current_v [[buffer(3)]],
                           const device solver_scalar* current_w [[buffer(4)]],
                           const device solver_scalar* pressure [[buffer(5)]],
                           const device solver_scalar* previous_advection [[buffer(6)]],
                           constant SolverParams& params [[buffer(7)]],
                           uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }

  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);

  const solver_scalar west_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_u, static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k), params) +
       load_field(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar east_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       load_field(current_u, static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar south_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_v, static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k), params) +
       load_field(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar north_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_v, static_cast<int>(i) - 1, static_cast<int>(j) + 1, static_cast<int>(k), params) +
       load_field(current_v, static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k), params));
  const solver_scalar back_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_w, static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k), params) +
       load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar front_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_w, static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k) + 1, params) +
       load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1, params));

  const solver_scalar west_flux =
      west_velocity * reconstruct_face_x(current_u, static_cast<int>(i), j, k, west_velocity, params);
  const solver_scalar east_flux =
      east_velocity * reconstruct_face_x(current_u, static_cast<int>(i) + 1, j, k, east_velocity, params);
  const solver_scalar south_flux =
      south_velocity * reconstruct_face_y(current_u, i, static_cast<int>(j), k, south_velocity, params);
  const solver_scalar north_flux =
      north_velocity * reconstruct_face_y(current_u, i, static_cast<int>(j) + 1, k, north_velocity, params);
  const solver_scalar back_flux =
      back_velocity * reconstruct_face_z(current_u, i, j, static_cast<int>(k), back_velocity, params);
  const solver_scalar front_flux =
      front_velocity * reconstruct_face_z(current_u, i, j, static_cast<int>(k) + 1, front_velocity, params);

  const solver_scalar advection =
      (east_flux - west_flux) / params.dx +
      (north_flux - south_flux) / params.dy +
      (front_flux - back_flux) / params.dz;
  const solver_scalar diffusion =
      params.viscosity *
      (second_derivative_x(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       second_derivative_y(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       second_derivative_z(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar factorized =
      mixed_second_derivative(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 0u, 1u, params) +
      mixed_second_derivative(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 0u, 2u, params) +
      mixed_second_derivative(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 1u, 2u, params) -
      params.alpha * triple_second_derivative(
          current_u,
          static_cast<int>(i),
          static_cast<int>(j),
          static_cast<int>(k),
          0u,
          1u,
          2u,
          params);

  const solver_scalar current_value =
      load_field(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params);
  solver_scalar rhs =
      current_value -
      params.dt * pressure_gradient_u(pressure, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
      static_cast<solver_scalar>(0.5) * params.dt * diffusion +
      params.alpha * params.alpha * factorized;

  if(params.use_ab2 != 0u) {
    rhs -= static_cast<solver_scalar>(1.5) * params.dt * advection;
    rhs += static_cast<solver_scalar>(0.5) * params.dt * previous_advection[gid];
  } else {
    rhs -= params.dt * advection;
  }

  advection_out[gid] = advection;
  rhs_out[gid] = rhs;
}

kernel void assemble_rhs_v(device solver_scalar* advection_out [[buffer(0)]],
                           device solver_scalar* rhs_out [[buffer(1)]],
                           const device solver_scalar* current_u [[buffer(2)]],
                           const device solver_scalar* current_v [[buffer(3)]],
                           const device solver_scalar* current_w [[buffer(4)]],
                           const device solver_scalar* pressure [[buffer(5)]],
                           const device solver_scalar* previous_advection [[buffer(6)]],
                           constant SolverParams& params [[buffer(7)]],
                           uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }

  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);

  const solver_scalar west_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_u, static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k), params) +
       load_field(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar east_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_u, static_cast<int>(i) + 1, static_cast<int>(j) - 1, static_cast<int>(k), params) +
       load_field(current_u, static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar south_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_v, static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k), params) +
       load_field(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar north_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       load_field(current_v, static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k), params));
  const solver_scalar back_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_w, static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k), params) +
       load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar front_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_w, static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k) + 1, params) +
       load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1, params));

  const solver_scalar west_flux =
      west_velocity * reconstruct_face_x(current_v, static_cast<int>(i), j, k, west_velocity, params);
  const solver_scalar east_flux =
      east_velocity * reconstruct_face_x(current_v, static_cast<int>(i) + 1, j, k, east_velocity, params);
  const solver_scalar south_flux =
      south_velocity * reconstruct_face_y(current_v, i, static_cast<int>(j), k, south_velocity, params);
  const solver_scalar north_flux =
      north_velocity * reconstruct_face_y(current_v, i, static_cast<int>(j) + 1, k, north_velocity, params);
  const solver_scalar back_flux =
      back_velocity * reconstruct_face_z(current_v, i, j, static_cast<int>(k), back_velocity, params);
  const solver_scalar front_flux =
      front_velocity * reconstruct_face_z(current_v, i, j, static_cast<int>(k) + 1, front_velocity, params);

  const solver_scalar advection =
      (east_flux - west_flux) / params.dx +
      (north_flux - south_flux) / params.dy +
      (front_flux - back_flux) / params.dz;
  const solver_scalar diffusion =
      params.viscosity *
      (second_derivative_x(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       second_derivative_y(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       second_derivative_z(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar factorized =
      mixed_second_derivative(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 0u, 1u, params) +
      mixed_second_derivative(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 0u, 2u, params) +
      mixed_second_derivative(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 1u, 2u, params) -
      params.alpha * triple_second_derivative(
          current_v,
          static_cast<int>(i),
          static_cast<int>(j),
          static_cast<int>(k),
          0u,
          1u,
          2u,
          params);

  const solver_scalar current_value =
      load_field(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params);
  solver_scalar rhs =
      current_value -
      params.dt * pressure_gradient_v(pressure, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
      static_cast<solver_scalar>(0.5) * params.dt * diffusion +
      params.alpha * params.alpha * factorized;

  if(params.use_ab2 != 0u) {
    rhs -= static_cast<solver_scalar>(1.5) * params.dt * advection;
    rhs += static_cast<solver_scalar>(0.5) * params.dt * previous_advection[gid];
  } else {
    rhs -= params.dt * advection;
  }

  advection_out[gid] = advection;
  rhs_out[gid] = rhs;
}

kernel void assemble_rhs_w(device solver_scalar* advection_out [[buffer(0)]],
                           device solver_scalar* rhs_out [[buffer(1)]],
                           const device solver_scalar* current_u [[buffer(2)]],
                           const device solver_scalar* current_v [[buffer(3)]],
                           const device solver_scalar* current_w [[buffer(4)]],
                           const device solver_scalar* pressure [[buffer(5)]],
                           const device solver_scalar* previous_advection [[buffer(6)]],
                           constant SolverParams& params [[buffer(7)]],
                           uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }

  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);

  const solver_scalar west_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1, params) +
       load_field(current_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar east_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_u, static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k) - 1, params) +
       load_field(current_u, static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar south_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1, params) +
       load_field(current_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar north_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_v, static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k) - 1, params) +
       load_field(current_v, static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k), params));
  const solver_scalar back_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1, params) +
       load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar front_velocity =
      static_cast<solver_scalar>(0.5) *
      (load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1, params));

  const solver_scalar west_flux =
      west_velocity * reconstruct_face_x(current_w, static_cast<int>(i), j, k, west_velocity, params);
  const solver_scalar east_flux =
      east_velocity * reconstruct_face_x(current_w, static_cast<int>(i) + 1, j, k, east_velocity, params);
  const solver_scalar south_flux =
      south_velocity * reconstruct_face_y(current_w, i, static_cast<int>(j), k, south_velocity, params);
  const solver_scalar north_flux =
      north_velocity * reconstruct_face_y(current_w, i, static_cast<int>(j) + 1, k, north_velocity, params);
  const solver_scalar back_flux =
      back_velocity * reconstruct_face_z(current_w, i, j, static_cast<int>(k), back_velocity, params);
  const solver_scalar front_flux =
      front_velocity * reconstruct_face_z(current_w, i, j, static_cast<int>(k) + 1, front_velocity, params);

  const solver_scalar advection =
      (east_flux - west_flux) / params.dx +
      (north_flux - south_flux) / params.dy +
      (front_flux - back_flux) / params.dz;
  const solver_scalar diffusion =
      params.viscosity *
      (second_derivative_x(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       second_derivative_y(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
       second_derivative_z(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params));
  const solver_scalar factorized =
      mixed_second_derivative(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 0u, 1u, params) +
      mixed_second_derivative(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 0u, 2u, params) +
      mixed_second_derivative(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), 1u, 2u, params) -
      params.alpha * triple_second_derivative(
          current_w,
          static_cast<int>(i),
          static_cast<int>(j),
          static_cast<int>(k),
          0u,
          1u,
          2u,
          params);

  const solver_scalar current_value =
      load_field(current_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params);
  solver_scalar rhs =
      current_value -
      params.dt * pressure_gradient_w(pressure, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
      static_cast<solver_scalar>(0.5) * params.dt * diffusion +
      params.alpha * params.alpha * factorized;

  if(params.use_ab2 != 0u) {
    rhs -= static_cast<solver_scalar>(1.5) * params.dt * advection;
    rhs += static_cast<solver_scalar>(0.5) * params.dt * previous_advection[gid];
  } else {
    rhs -= params.dt * advection;
  }

  advection_out[gid] = advection;
  rhs_out[gid] = rhs;
}

kernel void apply_factorized_operator(const device solver_scalar* input [[buffer(0)]],
                                      device solver_scalar* output [[buffer(1)]],
                                      constant SolverParams& params [[buffer(2)]],
                                      uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }

  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);
  output[gid] = factorized_operator_value(input, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params);
}

kernel void init_residual_and_direction(const device solver_scalar* rhs [[buffer(0)]],
                                        const device solver_scalar* applied [[buffer(1)]],
                                        device solver_scalar* residual [[buffer(2)]],
                                        device solver_scalar* direction [[buffer(3)]],
                                        constant uint& count [[buffer(4)]],
                                        uint gid [[thread_position_in_grid]]) {
  if(gid >= count) {
    return;
  }
  const solver_scalar value = rhs[gid] - applied[gid];
  residual[gid] = value;
  direction[gid] = value;
}

kernel void cg_update_solution_residual(device solver_scalar* solution [[buffer(0)]],
                                        device solver_scalar* residual [[buffer(1)]],
                                        const device solver_scalar* direction [[buffer(2)]],
                                        const device solver_scalar* applied [[buffer(3)]],
                                        constant solver_scalar& alpha [[buffer(4)]],
                                        constant uint& count [[buffer(5)]],
                                        uint gid [[thread_position_in_grid]]) {
  if(gid >= count) {
    return;
  }

  solution[gid] += alpha * direction[gid];
  residual[gid] -= alpha * applied[gid];
}

kernel void cg_update_direction(device solver_scalar* direction [[buffer(0)]],
                                const device solver_scalar* residual [[buffer(1)]],
                                constant solver_scalar& beta [[buffer(2)]],
                                constant uint& count [[buffer(3)]],
                                uint gid [[thread_position_in_grid]]) {
  if(gid >= count) {
    return;
  }
  direction[gid] = residual[gid] + beta * direction[gid];
}

kernel void subtract_scalar(device solver_scalar* values [[buffer(0)]],
                            constant solver_scalar& amount [[buffer(1)]],
                            constant uint& count [[buffer(2)]],
                            uint gid [[thread_position_in_grid]]) {
  if(gid >= count) {
    return;
  }
  values[gid] -= amount;
}

kernel void build_pressure_rhs(const device solver_scalar* predicted_u [[buffer(0)]],
                               const device solver_scalar* predicted_v [[buffer(1)]],
                               const device solver_scalar* predicted_w [[buffer(2)]],
                               device solver_scalar* rhs [[buffer(3)]],
                               constant SolverParams& params [[buffer(4)]],
                               uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }

  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);
  const solver_scalar divergence =
      (load_field(predicted_u, static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k), params) -
       load_field(predicted_u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params)) / params.dx +
      (load_field(predicted_v, static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k), params) -
       load_field(predicted_v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params)) / params.dy +
      (load_field(predicted_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1, params) -
       load_field(predicted_w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params)) / params.dz;
  rhs[gid] = -divergence / params.dt;
}

kernel void apply_poisson_operator(const device solver_scalar* input [[buffer(0)]],
                                   device solver_scalar* output [[buffer(1)]],
                                   constant SolverParams& params [[buffer(2)]],
                                   uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }

  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);
  const solver_scalar center = load_field(input, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params);
  const solver_scalar left = load_field(input, static_cast<int>(i) - 1, static_cast<int>(j), static_cast<int>(k), params);
  const solver_scalar right = load_field(input, static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k), params);
  const solver_scalar lower = load_field(input, static_cast<int>(i), static_cast<int>(j) - 1, static_cast<int>(k), params);
  const solver_scalar upper = load_field(input, static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k), params);
  const solver_scalar back = load_field(input, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) - 1, params);
  const solver_scalar front = load_field(input, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1, params);
  output[gid] =
      -(right - static_cast<solver_scalar>(2.0) * center + left) / (params.dx * params.dx) -
      (upper - static_cast<solver_scalar>(2.0) * center + lower) / (params.dy * params.dy) -
      (front - static_cast<solver_scalar>(2.0) * center + back) / (params.dz * params.dz);
}

kernel void correct_velocity_u(const device solver_scalar* predicted [[buffer(0)]],
                               const device solver_scalar* pressure_correction [[buffer(1)]],
                               device solver_scalar* corrected [[buffer(2)]],
                               constant SolverParams& params [[buffer(3)]],
                               uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }
  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);
  corrected[gid] =
      predicted[gid] -
      params.dt * pressure_gradient_u(
          pressure_correction,
          static_cast<int>(i),
          static_cast<int>(j),
          static_cast<int>(k),
          params);
}

kernel void correct_velocity_v(const device solver_scalar* predicted [[buffer(0)]],
                               const device solver_scalar* pressure_correction [[buffer(1)]],
                               device solver_scalar* corrected [[buffer(2)]],
                               constant SolverParams& params [[buffer(3)]],
                               uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }
  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);
  corrected[gid] =
      predicted[gid] -
      params.dt * pressure_gradient_v(
          pressure_correction,
          static_cast<int>(i),
          static_cast<int>(j),
          static_cast<int>(k),
          params);
}

kernel void correct_velocity_w(const device solver_scalar* predicted [[buffer(0)]],
                               const device solver_scalar* pressure_correction [[buffer(1)]],
                               device solver_scalar* corrected [[buffer(2)]],
                               constant SolverParams& params [[buffer(3)]],
                               uint gid [[thread_position_in_grid]]) {
  const uint count = params.nx * params.ny * params.nz;
  if(gid >= count) {
    return;
  }
  uint i = 0;
  uint j = 0;
  uint k = 0;
  unpack_index(gid, params, i, j, k);
  corrected[gid] =
      predicted[gid] -
      params.dt * pressure_gradient_w(
          pressure_correction,
          static_cast<int>(i),
          static_cast<int>(j),
          static_cast<int>(k),
          params);
}

kernel void add_pressure_correction(device solver_scalar* pressure_total [[buffer(0)]],
                                    const device solver_scalar* pressure_correction [[buffer(1)]],
                                    constant uint& count [[buffer(2)]],
                                    uint gid [[thread_position_in_grid]]) {
  if(gid >= count) {
    return;
  }
  pressure_total[gid] += pressure_correction[gid];
}

kernel void dot_partial(const device solver_scalar* left [[buffer(0)]],
                        const device solver_scalar* right [[buffer(1)]],
                        device solver_scalar* partial [[buffer(2)]],
                        constant uint& count [[buffer(3)]],
                        uint gid [[thread_position_in_grid]],
                        uint lid [[thread_index_in_threadgroup]],
                        uint group_id [[threadgroup_position_in_grid]]) {
  threadgroup solver_scalar scratch[kReductionWidth];
  scratch[lid] = gid < count ? left[gid] * right[gid] : static_cast<solver_scalar>(0.0);
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for(uint stride = kReductionWidth / 2u; stride > 0u; stride >>= 1u) {
    if(lid < stride) {
      scratch[lid] += scratch[lid + stride];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }

  if(lid == 0u) {
    partial[group_id] = scratch[0];
  }
}

kernel void sum_partial(const device solver_scalar* input [[buffer(0)]],
                        device solver_scalar* partial [[buffer(1)]],
                        constant uint& count [[buffer(2)]],
                        uint gid [[thread_position_in_grid]],
                        uint lid [[thread_index_in_threadgroup]],
                        uint group_id [[threadgroup_position_in_grid]]) {
  threadgroup solver_scalar scratch[kReductionWidth];
  scratch[lid] = gid < count ? input[gid] : static_cast<solver_scalar>(0.0);
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for(uint stride = kReductionWidth / 2u; stride > 0u; stride >>= 1u) {
    if(lid < stride) {
      scratch[lid] += scratch[lid + stride];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }

  if(lid == 0u) {
    partial[group_id] = scratch[0];
  }
}

kernel void max_abs_diff_partial(const device solver_scalar* left [[buffer(0)]],
                                 const device solver_scalar* right [[buffer(1)]],
                                 device solver_scalar* partial [[buffer(2)]],
                                 constant uint& count [[buffer(3)]],
                                 uint gid [[thread_position_in_grid]],
                                 uint lid [[thread_index_in_threadgroup]],
                                 uint group_id [[threadgroup_position_in_grid]]) {
  threadgroup solver_scalar scratch[kReductionWidth];
  scratch[lid] = gid < count ? fabs(left[gid] - right[gid]) : static_cast<solver_scalar>(0.0);
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for(uint stride = kReductionWidth / 2u; stride > 0u; stride >>= 1u) {
    if(lid < stride) {
      scratch[lid] = max(scratch[lid], scratch[lid + stride]);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }

  if(lid == 0u) {
    partial[group_id] = scratch[0];
  }
}

kernel void cfl_partial(const device solver_scalar* u [[buffer(0)]],
                        const device solver_scalar* v [[buffer(1)]],
                        const device solver_scalar* w [[buffer(2)]],
                        device solver_scalar* partial [[buffer(3)]],
                        constant SolverParams& params [[buffer(4)]],
                        uint gid [[thread_position_in_grid]],
                        uint lid [[thread_index_in_threadgroup]],
                        uint group_id [[threadgroup_position_in_grid]]) {
  threadgroup solver_scalar scratch[kReductionWidth];
  const uint count = params.nx * params.ny * params.nz;
  solver_scalar value = static_cast<solver_scalar>(0.0);

  if(gid < count) {
    uint i = 0;
    uint j = 0;
    uint k = 0;
    unpack_index(gid, params, i, j, k);
    const solver_scalar u_center =
        static_cast<solver_scalar>(0.5) *
        (load_field(u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
         load_field(u, static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k), params));
    const solver_scalar v_center =
        static_cast<solver_scalar>(0.5) *
        (load_field(v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
         load_field(v, static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k), params));
    const solver_scalar w_center =
        static_cast<solver_scalar>(0.5) *
        (load_field(w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params) +
         load_field(w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1, params));
    value = params.dt *
            (fabs(u_center) / params.dx + fabs(v_center) / params.dy + fabs(w_center) / params.dz);
  }

  scratch[lid] = value;
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for(uint stride = kReductionWidth / 2u; stride > 0u; stride >>= 1u) {
    if(lid < stride) {
      scratch[lid] = max(scratch[lid], scratch[lid + stride]);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }

  if(lid == 0u) {
    partial[group_id] = scratch[0];
  }
}

kernel void divergence_square_partial(const device solver_scalar* u [[buffer(0)]],
                                      const device solver_scalar* v [[buffer(1)]],
                                      const device solver_scalar* w [[buffer(2)]],
                                      device solver_scalar* partial [[buffer(3)]],
                                      constant SolverParams& params [[buffer(4)]],
                                      uint gid [[thread_position_in_grid]],
                                      uint lid [[thread_index_in_threadgroup]],
                                      uint group_id [[threadgroup_position_in_grid]]) {
  threadgroup solver_scalar scratch[kReductionWidth];
  const uint count = params.nx * params.ny * params.nz;
  solver_scalar value = static_cast<solver_scalar>(0.0);

  if(gid < count) {
    uint i = 0;
    uint j = 0;
    uint k = 0;
    unpack_index(gid, params, i, j, k);
    const solver_scalar divergence =
        (load_field(u, static_cast<int>(i) + 1, static_cast<int>(j), static_cast<int>(k), params) -
         load_field(u, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params)) / params.dx +
        (load_field(v, static_cast<int>(i), static_cast<int>(j) + 1, static_cast<int>(k), params) -
         load_field(v, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params)) / params.dy +
        (load_field(w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k) + 1, params) -
         load_field(w, static_cast<int>(i), static_cast<int>(j), static_cast<int>(k), params)) / params.dz;
    value = divergence * divergence;
  }

  scratch[lid] = value;
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for(uint stride = kReductionWidth / 2u; stride > 0u; stride >>= 1u) {
    if(lid < stride) {
      scratch[lid] += scratch[lid + stride];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }

  if(lid == 0u) {
    partial[group_id] = scratch[0];
  }
}
