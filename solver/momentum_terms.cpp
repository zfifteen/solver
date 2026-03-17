#include "solver/momentum_terms.hpp"

#include "operators/discrete_operators.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace solver {

namespace {

constexpr double kSmallDenominator = 1.0e-14;

void require_same_velocity_layouts(const VelocityField& left,
                                   const VelocityField& right,
                                   const char* operation_name) {
  if(!left.x.layout().same_shape_as(right.x.layout()) ||
     !left.y.layout().same_shape_as(right.y.layout()) ||
     !left.z.layout().same_shape_as(right.z.layout())) {
    throw std::invalid_argument(std::string(operation_name) + ": incompatible velocity layouts");
  }
}

double limiter_value(const FluxLimiter limiter, const double ratio) {
  switch(limiter) {
    case FluxLimiter::van_leer:
      return (ratio + std::abs(ratio)) / (1.0 + std::abs(ratio));
  }

  __builtin_unreachable();
}

double safe_ratio(const double numerator, const double denominator) {
  if(std::abs(denominator) < kSmallDenominator) {
    return 0.0;
  }
  return numerator / denominator;
}

double reconstruct_face_x(const StructuredField& field,
                          const int i_right,
                          const int j,
                          const int k,
                          const double face_velocity,
                          const AdvectionOptions& options) {
  const double left = field(i_right - 1, j, k);
  const double right = field(i_right, j, k);

  if(options.scheme == AdvectionScheme::central) {
    return 0.5 * (left + right);
  }

  if(options.scheme == AdvectionScheme::upwind || std::abs(face_velocity) < kSmallDenominator) {
    return face_velocity >= 0.0 ? left : right;
  }

  const int max_i = field.layout().storage_extent().nx - 1;
  if(face_velocity >= 0.0) {
    if(i_right - 2 < 0) {
      return left;
    }

    const double delta_down = right - left;
    const double ratio = safe_ratio(left - field(i_right - 2, j, k), delta_down);
    return left + 0.5 * limiter_value(options.limiter, ratio) * delta_down;
  }

  if(i_right + 1 > max_i) {
    return right;
  }

  const double delta = right - left;
  const double ratio = safe_ratio(field(i_right + 1, j, k) - right, delta);
  return right - 0.5 * limiter_value(options.limiter, ratio) * delta;
}

double reconstruct_face_y(const StructuredField& field,
                          const int i,
                          const int j_right,
                          const int k,
                          const double face_velocity,
                          const AdvectionOptions& options) {
  const double lower = field(i, j_right - 1, k);
  const double upper = field(i, j_right, k);

  if(options.scheme == AdvectionScheme::central) {
    return 0.5 * (lower + upper);
  }

  if(options.scheme == AdvectionScheme::upwind || std::abs(face_velocity) < kSmallDenominator) {
    return face_velocity >= 0.0 ? lower : upper;
  }

  const int max_j = field.layout().storage_extent().ny - 1;
  if(face_velocity >= 0.0) {
    if(j_right - 2 < 0) {
      return lower;
    }

    const double delta_down = upper - lower;
    const double ratio = safe_ratio(lower - field(i, j_right - 2, k), delta_down);
    return lower + 0.5 * limiter_value(options.limiter, ratio) * delta_down;
  }

  if(j_right + 1 > max_j) {
    return upper;
  }

  const double delta = upper - lower;
  const double ratio = safe_ratio(field(i, j_right + 1, k) - upper, delta);
  return upper - 0.5 * limiter_value(options.limiter, ratio) * delta;
}

double reconstruct_face_z(const StructuredField& field,
                          const int i,
                          const int j,
                          const int k_right,
                          const double face_velocity,
                          const AdvectionOptions& options) {
  const double back = field(i, j, k_right - 1);
  const double front = field(i, j, k_right);

  if(options.scheme == AdvectionScheme::central) {
    return 0.5 * (back + front);
  }

  if(options.scheme == AdvectionScheme::upwind || std::abs(face_velocity) < kSmallDenominator) {
    return face_velocity >= 0.0 ? back : front;
  }

  const int max_k = field.layout().storage_extent().nz - 1;
  if(face_velocity >= 0.0) {
    if(k_right - 2 < 0) {
      return back;
    }

    const double delta_down = front - back;
    const double ratio = safe_ratio(back - field(i, j, k_right - 2), delta_down);
    return back + 0.5 * limiter_value(options.limiter, ratio) * delta_down;
  }

  if(k_right + 1 > max_k) {
    return front;
  }

  const double delta = front - back;
  const double ratio = safe_ratio(field(i, j, k_right + 1) - front, delta);
  return front - 0.5 * limiter_value(options.limiter, ratio) * delta;
}

void compute_u_advection(const VelocityField& velocity,
                         const AdvectionOptions& options,
                         FaceField& output) {
  const IndexRange3D active = velocity.x.layout().active_range();
  const Grid& grid = velocity.x.layout().grid();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double west_velocity = 0.5 * (velocity.x(i - 1, j, k) + velocity.x(i, j, k));
        const double east_velocity = 0.5 * (velocity.x(i, j, k) + velocity.x(i + 1, j, k));
        const double south_velocity = 0.5 * (velocity.y(i - 1, j, k) + velocity.y(i, j, k));
        const double north_velocity = 0.5 * (velocity.y(i - 1, j + 1, k) + velocity.y(i, j + 1, k));
        const double back_velocity = 0.5 * (velocity.z(i - 1, j, k) + velocity.z(i, j, k));
        const double front_velocity =
            0.5 * (velocity.z(i - 1, j, k + 1) + velocity.z(i, j, k + 1));

        const double west_flux =
            west_velocity * reconstruct_face_x(velocity.x, i, j, k, west_velocity, options);
        const double east_flux =
            east_velocity * reconstruct_face_x(velocity.x, i + 1, j, k, east_velocity, options);
        const double south_flux =
            south_velocity * reconstruct_face_y(velocity.x, i, j, k, south_velocity, options);
        const double north_flux =
            north_velocity * reconstruct_face_y(velocity.x, i, j + 1, k, north_velocity, options);
        const double back_flux =
            back_velocity * reconstruct_face_z(velocity.x, i, j, k, back_velocity, options);
        const double front_flux =
            front_velocity * reconstruct_face_z(velocity.x, i, j, k + 1, front_velocity, options);

        output(i, j, k) = (east_flux - west_flux) / grid.dx +
                          (north_flux - south_flux) / grid.dy +
                          (front_flux - back_flux) / grid.dz;
      }
    }
  }
}

void compute_v_advection(const VelocityField& velocity,
                         const AdvectionOptions& options,
                         FaceField& output) {
  const IndexRange3D active = velocity.y.layout().active_range();
  const Grid& grid = velocity.y.layout().grid();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double west_velocity = 0.5 * (velocity.x(i, j - 1, k) + velocity.x(i, j, k));
        const double east_velocity = 0.5 * (velocity.x(i + 1, j - 1, k) + velocity.x(i + 1, j, k));
        const double south_velocity = 0.5 * (velocity.y(i, j - 1, k) + velocity.y(i, j, k));
        const double north_velocity = 0.5 * (velocity.y(i, j, k) + velocity.y(i, j + 1, k));
        const double back_velocity = 0.5 * (velocity.z(i, j - 1, k) + velocity.z(i, j, k));
        const double front_velocity =
            0.5 * (velocity.z(i, j - 1, k + 1) + velocity.z(i, j, k + 1));

        const double west_flux =
            west_velocity * reconstruct_face_x(velocity.y, i, j, k, west_velocity, options);
        const double east_flux =
            east_velocity * reconstruct_face_x(velocity.y, i + 1, j, k, east_velocity, options);
        const double south_flux =
            south_velocity * reconstruct_face_y(velocity.y, i, j, k, south_velocity, options);
        const double north_flux =
            north_velocity * reconstruct_face_y(velocity.y, i, j + 1, k, north_velocity, options);
        const double back_flux =
            back_velocity * reconstruct_face_z(velocity.y, i, j, k, back_velocity, options);
        const double front_flux =
            front_velocity * reconstruct_face_z(velocity.y, i, j, k + 1, front_velocity, options);

        output(i, j, k) = (east_flux - west_flux) / grid.dx +
                          (north_flux - south_flux) / grid.dy +
                          (front_flux - back_flux) / grid.dz;
      }
    }
  }
}

void compute_w_advection(const VelocityField& velocity,
                         const AdvectionOptions& options,
                         FaceField& output) {
  const IndexRange3D active = velocity.z.layout().active_range();
  const Grid& grid = velocity.z.layout().grid();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double west_velocity = 0.5 * (velocity.x(i, j, k - 1) + velocity.x(i, j, k));
        const double east_velocity =
            0.5 * (velocity.x(i + 1, j, k - 1) + velocity.x(i + 1, j, k));
        const double south_velocity = 0.5 * (velocity.y(i, j, k - 1) + velocity.y(i, j, k));
        const double north_velocity =
            0.5 * (velocity.y(i, j + 1, k - 1) + velocity.y(i, j + 1, k));
        const double back_velocity = 0.5 * (velocity.z(i, j, k - 1) + velocity.z(i, j, k));
        const double front_velocity = 0.5 * (velocity.z(i, j, k) + velocity.z(i, j, k + 1));

        const double west_flux =
            west_velocity * reconstruct_face_x(velocity.z, i, j, k, west_velocity, options);
        const double east_flux =
            east_velocity * reconstruct_face_x(velocity.z, i + 1, j, k, east_velocity, options);
        const double south_flux =
            south_velocity * reconstruct_face_y(velocity.z, i, j, k, south_velocity, options);
        const double north_flux =
            north_velocity * reconstruct_face_y(velocity.z, i, j + 1, k, north_velocity, options);
        const double back_flux =
            back_velocity * reconstruct_face_z(velocity.z, i, j, k, back_velocity, options);
        const double front_flux =
            front_velocity * reconstruct_face_z(velocity.z, i, j, k + 1, front_velocity, options);

        output(i, j, k) = (east_flux - west_flux) / grid.dx +
                          (north_flux - south_flux) / grid.dy +
                          (front_flux - back_flux) / grid.dz;
      }
    }
  }
}

void scale_active(StructuredField& field, const double scale) {
  const IndexRange3D active = field.layout().active_range();
  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        field(i, j, k) *= scale;
      }
    }
  }
}

}  // namespace

std::string to_string(const AdvectionScheme scheme) {
  switch(scheme) {
    case AdvectionScheme::tvd:
      return "tvd";
    case AdvectionScheme::upwind:
      return "upwind";
    case AdvectionScheme::central:
      return "central";
  }

  __builtin_unreachable();
}

std::string to_string(const FluxLimiter limiter) {
  switch(limiter) {
    case FluxLimiter::van_leer:
      return "van_leer";
  }

  __builtin_unreachable();
}

std::string describe(const AdvectionOptions& options) {
  return "scheme=" + to_string(options.scheme) + ", limiter=" + to_string(options.limiter);
}

void compute_advection_term(const VelocityField& velocity,
                            const AdvectionOptions& options,
                            VelocityField& advection) {
  require_same_velocity_layouts(velocity, advection, "compute_advection_term");
  compute_u_advection(velocity, options, advection.x);
  compute_v_advection(velocity, options, advection.y);
  compute_w_advection(velocity, options, advection.z);
}

void compute_diffusion_term(const VelocityField& velocity,
                            const double viscosity,
                            VelocityField& diffusion) {
  require_same_velocity_layouts(velocity, diffusion, "compute_diffusion_term");

  operators::compute_laplacian(velocity, diffusion);
  scale_active(diffusion.x, viscosity);
  scale_active(diffusion.y, viscosity);
  scale_active(diffusion.z, viscosity);
}

CflDiagnostics compute_advective_cfl(const VelocityField& velocity, const double dt) {
  if(dt < 0.0) {
    throw std::invalid_argument("compute_advective_cfl expects a non-negative dt");
  }

  const Grid& grid = velocity.x.layout().grid();
  const IndexRange3D cells = FieldLayout::cell_centered(grid).active_range();
  CflDiagnostics diagnostics{.dt = dt};

  for(int k = cells.k_begin; k < cells.k_end; ++k) {
    for(int j = cells.j_begin; j < cells.j_end; ++j) {
      for(int i = cells.i_begin; i < cells.i_end; ++i) {
        const double u_center = 0.5 * (velocity.x(i, j, k) + velocity.x(i + 1, j, k));
        const double v_center = 0.5 * (velocity.y(i, j, k) + velocity.y(i, j + 1, k));
        const double w_center = 0.5 * (velocity.z(i, j, k) + velocity.z(i, j, k + 1));

        diagnostics.max_u = std::max(diagnostics.max_u, std::abs(u_center));
        diagnostics.max_v = std::max(diagnostics.max_v, std::abs(v_center));
        diagnostics.max_w = std::max(diagnostics.max_w, std::abs(w_center));

        const double local_cfl = dt * (std::abs(u_center) / grid.dx + std::abs(v_center) / grid.dy +
                                       std::abs(w_center) / grid.dz);
        diagnostics.max_cfl = std::max(diagnostics.max_cfl, local_cfl);
      }
    }
  }

  return diagnostics;
}

}  // namespace solver
