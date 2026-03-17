#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "metal/taylor_green_backend.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace solver::metal {

namespace {

#if SOLVER_METAL_USE_FLOAT
using MetalScalar = float;
constexpr double kCgTolerance = 1.0e-6;
#else
using MetalScalar = double;
constexpr double kCgTolerance = 1.0e-12;
#endif

constexpr NSUInteger kThreadsPerThreadgroup = 256;

struct SolverParams {
  uint32_t nx = 0;
  uint32_t ny = 0;
  uint32_t nz = 0;
  uint32_t use_ab2 = 0;
  uint32_t advection_scheme = 0;
  MetalScalar dx = 0;
  MetalScalar dy = 0;
  MetalScalar dz = 0;
  MetalScalar viscosity = 0;
  MetalScalar dt = 0;
  MetalScalar alpha = 0;
};

static_assert(std::is_standard_layout_v<SolverParams>);

std::string ns_string_to_std(NSString* value) {
  return value == nil ? std::string{} : std::string([value UTF8String]);
}

std::string build_error_message(const std::string& prefix, NSError* error) {
  std::ostringstream builder;
  builder << prefix;
  if(error != nil) {
    builder << ": " << ns_string_to_std(error.localizedDescription);
  }
  return builder.str();
}

void require_or_throw(const bool condition, const std::string& message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

Grid make_grid_from_config(const TaylorGreenConfig& config) {
  const double domain_length = 2.0 * std::acos(-1.0);
  return Grid{
      config.nx,
      config.ny,
      config.nz,
      domain_length / static_cast<double>(config.nx),
      domain_length / static_cast<double>(config.ny),
      domain_length / static_cast<double>(config.nz),
      1,
  };
}

double compute_dt(const TaylorGreenConfig& config) {
  const Grid grid = make_grid_from_config(config);
  const double inverse_sum = 1.0 / grid.dx + 1.0 / grid.dy + 1.0 / grid.dz;
  const double dt_upper = config.cfl_limit / inverse_sum;
  const int step_count =
      std::max(1, static_cast<int>(std::ceil(config.final_time / dt_upper)));
  return config.final_time / static_cast<double>(step_count);
}

uint32_t partial_count_for(const uint32_t count) {
  return (count + static_cast<uint32_t>(kThreadsPerThreadgroup) - 1u) /
         static_cast<uint32_t>(kThreadsPerThreadgroup);
}

template <typename T>
T* buffer_data(id<MTLBuffer> buffer) {
  return static_cast<T*>([buffer contents]);
}

class MetalRunner {
 public:
  explicit MetalRunner(const TaylorGreenConfig& config)
      : config_(config),
        grid_(make_grid_from_config(config_)),
        step_count_(static_cast<int>(std::ceil(config_.final_time / compute_dt(config_)))),
        dt_(compute_dt(config_)),
        params_(make_params()),
        face_count_(static_cast<uint32_t>(config_.nx * config_.ny * config_.nz)),
        cell_count_(face_count_),
        partial_count_(partial_count_for(cell_count_)) {
    if(config_.backend != ExecutionBackend::metal) {
      throw std::invalid_argument("Metal runner expects backend=metal");
    }
    if(config_.nz <= 1) {
      throw std::invalid_argument("Metal Taylor-Green backend supports 3D periodic cases only");
    }

    @autoreleasepool {
      device_ = MTLCreateSystemDefaultDevice();
      require_or_throw(device_ != nil, "Metal backend could not find a default system device");
      command_queue_ = [device_ newCommandQueue];
      require_or_throw(command_queue_ != nil, "Metal backend could not create a command queue");

      NSError* library_error = nil;
      NSString* library_path = [NSString stringWithUTF8String:SOLVER_METAL_LIBRARY_PATH];
      NSURL* library_url = [NSURL fileURLWithPath:library_path];
      library_ = [device_ newLibraryWithURL:library_url error:&library_error];
      require_or_throw(library_ != nil,
                       build_error_message("Metal backend could not load the generated metallib",
                                           library_error));

      load_pipelines();
      allocate_buffers();
      zero_buffer(previous_advection_[0]);
      zero_buffer(previous_advection_[1]);
      zero_buffer(previous_advection_[2]);
      zero_buffer(pressure_correction_);
      zero_buffer(pressure_residual_);
      zero_buffer(pressure_direction_);
      zero_buffer(pressure_applied_);
      zero_buffer(face_rhs_);
      zero_buffer(face_residual_);
      zero_buffer(face_direction_);
      zero_buffer(face_applied_);
    }
  }

  TaylorGreenMetalRun run() {
    const auto started = std::chrono::steady_clock::now();
    initialize_state();

    TaylorGreenStepMetrics metrics{};
    metrics.dt = dt_;
    bool has_previous_advection = false;
    double max_divergence_seen = 0.0;

    for(int step = 0; step < step_count_; ++step) {
      params_.use_ab2 = has_previous_advection ? 1u : 0u;

      assemble_component(0u);
      solve_face_system(next_velocity_[0]);
      assemble_component(1u);
      solve_face_system(next_velocity_[1]);
      assemble_component(2u);
      solve_face_system(next_velocity_[2]);

      build_pressure_rhs();
      const MetalScalar rhs_mean = reduce_sum(pressure_rhs_) / static_cast<MetalScalar>(cell_count_);
      if(std::abs(static_cast<double>(rhs_mean)) > 0.0) {
        subtract_scalar_from_buffer(pressure_rhs_, rhs_mean, cell_count_);
      }

      const auto pressure_stats = solve_pressure_system();

      correct_velocity(0u);
      correct_velocity(1u);
      correct_velocity(2u);
      add_pressure_correction();

      const double divergence_l2 = compute_divergence_l2(next_velocity_);
      const double max_cfl = compute_cfl(next_velocity_);
      const double max_change = std::max(
          {compute_max_abs_difference(current_velocity_[0], next_velocity_[0]),
           compute_max_abs_difference(current_velocity_[1], next_velocity_[1]),
           compute_max_abs_difference(current_velocity_[2], next_velocity_[2])});

      max_divergence_seen = std::max(max_divergence_seen, divergence_l2);
      metrics = TaylorGreenStepMetrics{
          .step = step + 1,
          .time = static_cast<double>(step + 1) * dt_,
          .dt = dt_,
          .max_cfl = max_cfl,
          .max_velocity_change = max_change,
          .divergence_l2 = divergence_l2,
          .max_divergence_l2 = max_divergence_seen,
          .pressure_iterations = pressure_stats.iterations,
          .pressure_relative_residual = pressure_stats.relative_residual,
      };

      if(!std::isfinite(divergence_l2) || !std::isfinite(max_change) || !std::isfinite(max_cfl)) {
        throw std::runtime_error("Metal Taylor-Green backend produced a non-finite state");
      }

      std::swap(current_velocity_, next_velocity_);
      std::swap(previous_advection_, current_advection_);
      has_previous_advection = true;
    }

    TaylorGreenMetalRun run{
        .state = unpack_state(metrics, has_previous_advection),
        .device_name = ns_string_to_std([device_ name]),
        .elapsed_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count(),
    };
    return run;
  }

 private:
  struct PressureSolveStats {
    int iterations = 0;
    double relative_residual = 0.0;
  };

  SolverParams make_params() const {
    return SolverParams{
        .nx = static_cast<uint32_t>(config_.nx),
        .ny = static_cast<uint32_t>(config_.ny),
        .nz = static_cast<uint32_t>(config_.nz),
        .use_ab2 = 0u,
        .advection_scheme = static_cast<uint32_t>(config_.advection.scheme),
        .dx = static_cast<MetalScalar>(grid_.dx),
        .dy = static_cast<MetalScalar>(grid_.dy),
        .dz = static_cast<MetalScalar>(grid_.dz),
        .viscosity = static_cast<MetalScalar>(config_.viscosity),
        .dt = static_cast<MetalScalar>(dt_),
        .alpha = static_cast<MetalScalar>(0.5 * config_.viscosity * dt_),
    };
  }

  id<MTLComputePipelineState> make_pipeline(NSString* name) {
    NSError* error = nil;
    id<MTLFunction> function = [library_ newFunctionWithName:name];
    require_or_throw(function != nil,
                     "Metal backend could not find kernel '" + ns_string_to_std(name) + "'");
    id<MTLComputePipelineState> pipeline =
        [device_ newComputePipelineStateWithFunction:function error:&error];
    require_or_throw(
        pipeline != nil,
        build_error_message("Metal backend could not compile pipeline '" + ns_string_to_std(name) + "'",
                            error));
    return pipeline;
  }

  void load_pipelines() {
    init_pipeline_ = make_pipeline(@"init_taylor_green_3d");
    assemble_u_pipeline_ = make_pipeline(@"assemble_rhs_u");
    assemble_v_pipeline_ = make_pipeline(@"assemble_rhs_v");
    assemble_w_pipeline_ = make_pipeline(@"assemble_rhs_w");
    factorized_operator_pipeline_ = make_pipeline(@"apply_factorized_operator");
    init_residual_pipeline_ = make_pipeline(@"init_residual_and_direction");
    update_solution_pipeline_ = make_pipeline(@"cg_update_solution_residual");
    update_direction_pipeline_ = make_pipeline(@"cg_update_direction");
    subtract_scalar_pipeline_ = make_pipeline(@"subtract_scalar");
    pressure_rhs_pipeline_ = make_pipeline(@"build_pressure_rhs");
    poisson_operator_pipeline_ = make_pipeline(@"apply_poisson_operator");
    correct_u_pipeline_ = make_pipeline(@"correct_velocity_u");
    correct_v_pipeline_ = make_pipeline(@"correct_velocity_v");
    correct_w_pipeline_ = make_pipeline(@"correct_velocity_w");
    add_pressure_pipeline_ = make_pipeline(@"add_pressure_correction");
    dot_partial_pipeline_ = make_pipeline(@"dot_partial");
    sum_partial_pipeline_ = make_pipeline(@"sum_partial");
    max_abs_diff_pipeline_ = make_pipeline(@"max_abs_diff_partial");
    cfl_partial_pipeline_ = make_pipeline(@"cfl_partial");
    divergence_partial_pipeline_ = make_pipeline(@"divergence_square_partial");
  }

  id<MTLBuffer> make_buffer(const NSUInteger bytes) {
    id<MTLBuffer> buffer =
        [device_ newBufferWithLength:bytes options:MTLResourceStorageModeShared];
    require_or_throw(buffer != nil, "Metal backend failed to allocate a shared buffer");
    return buffer;
  }

  id<MTLBuffer> make_scalar_buffer(const uint32_t count) {
    return make_buffer(static_cast<NSUInteger>(count) * sizeof(MetalScalar));
  }

  void allocate_buffers() {
    for(std::size_t component = 0; component < 3; ++component) {
      current_velocity_[component] = make_scalar_buffer(face_count_);
      next_velocity_[component] = make_scalar_buffer(face_count_);
      previous_advection_[component] = make_scalar_buffer(face_count_);
      current_advection_[component] = make_scalar_buffer(face_count_);
    }

    pressure_total_ = make_scalar_buffer(cell_count_);
    pressure_rhs_ = make_scalar_buffer(cell_count_);
    pressure_correction_ = make_scalar_buffer(cell_count_);
    pressure_residual_ = make_scalar_buffer(cell_count_);
    pressure_direction_ = make_scalar_buffer(cell_count_);
    pressure_applied_ = make_scalar_buffer(cell_count_);

    face_rhs_ = make_scalar_buffer(face_count_);
    face_residual_ = make_scalar_buffer(face_count_);
    face_direction_ = make_scalar_buffer(face_count_);
    face_applied_ = make_scalar_buffer(face_count_);

    partial_buffer_ = make_scalar_buffer(partial_count_);
  }

  template <typename Encode>
  void dispatch_1d(id<MTLComputePipelineState> pipeline,
                   const uint32_t count,
                   Encode&& encode) {
    @autoreleasepool {
      id<MTLCommandBuffer> command_buffer = [command_queue_ commandBuffer];
      id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
      [encoder setComputePipelineState:pipeline];
      encode(encoder);
      const MTLSize threads_per_group = MTLSizeMake(kThreadsPerThreadgroup, 1, 1);
      const MTLSize threads_per_grid = MTLSizeMake(count, 1, 1);
      [encoder dispatchThreads:threads_per_grid threadsPerThreadgroup:threads_per_group];
      [encoder endEncoding];
      [command_buffer commit];
      [command_buffer waitUntilCompleted];
      require_or_throw(command_buffer.error == nil,
                       build_error_message("Metal compute dispatch failed", command_buffer.error));
    }
  }

  void zero_buffer(id<MTLBuffer> buffer) {
    @autoreleasepool {
      id<MTLCommandBuffer> command_buffer = [command_queue_ commandBuffer];
      id<MTLBlitCommandEncoder> encoder = [command_buffer blitCommandEncoder];
      [encoder fillBuffer:buffer range:NSMakeRange(0, [buffer length]) value:0];
      [encoder endEncoding];
      [command_buffer commit];
      [command_buffer waitUntilCompleted];
      require_or_throw(command_buffer.error == nil,
                       build_error_message("Metal blit fill failed", command_buffer.error));
    }
  }

  void copy_buffer(id<MTLBuffer> source, id<MTLBuffer> destination, const NSUInteger bytes) {
    @autoreleasepool {
      id<MTLCommandBuffer> command_buffer = [command_queue_ commandBuffer];
      id<MTLBlitCommandEncoder> encoder = [command_buffer blitCommandEncoder];
      [encoder copyFromBuffer:source
                 sourceOffset:0
                     toBuffer:destination
            destinationOffset:0
                         size:bytes];
      [encoder endEncoding];
      [command_buffer commit];
      [command_buffer waitUntilCompleted];
      require_or_throw(command_buffer.error == nil,
                       build_error_message("Metal blit copy failed", command_buffer.error));
    }
  }

  void initialize_state() {
    dispatch_1d(
        init_pipeline_,
        cell_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:current_velocity_[0] offset:0 atIndex:0];
          [encoder setBuffer:current_velocity_[1] offset:0 atIndex:1];
          [encoder setBuffer:current_velocity_[2] offset:0 atIndex:2];
          [encoder setBuffer:pressure_total_ offset:0 atIndex:3];
          [encoder setBytes:&params_ length:sizeof(params_) atIndex:4];
        });
  }

  void assemble_component(const uint32_t component) {
    id<MTLComputePipelineState> pipeline = assemble_u_pipeline_;
    switch(component) {
      case 0u:
        pipeline = assemble_u_pipeline_;
        break;
      case 1u:
        pipeline = assemble_v_pipeline_;
        break;
      case 2u:
        pipeline = assemble_w_pipeline_;
        break;
    }

    dispatch_1d(
        pipeline,
        face_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:current_advection_[component] offset:0 atIndex:0];
          [encoder setBuffer:face_rhs_ offset:0 atIndex:1];
          [encoder setBuffer:current_velocity_[0] offset:0 atIndex:2];
          [encoder setBuffer:current_velocity_[1] offset:0 atIndex:3];
          [encoder setBuffer:current_velocity_[2] offset:0 atIndex:4];
          [encoder setBuffer:pressure_total_ offset:0 atIndex:5];
          [encoder setBuffer:previous_advection_[component] offset:0 atIndex:6];
          [encoder setBytes:&params_ length:sizeof(params_) atIndex:7];
        });
  }

  MetalScalar reduce_sum(id<MTLBuffer> buffer) {
    const uint32_t count = static_cast<uint32_t>([buffer length] / sizeof(MetalScalar));
    const uint32_t active_count =
        std::min(count, std::max(cell_count_, face_count_));
    dispatch_1d(
        sum_partial_pipeline_,
        active_count,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:buffer offset:0 atIndex:0];
          [encoder setBuffer:partial_buffer_ offset:0 atIndex:1];
          [encoder setBytes:&active_count length:sizeof(active_count) atIndex:2];
        });

    const MetalScalar* partial = buffer_data<MetalScalar>(partial_buffer_);
    MetalScalar sum = static_cast<MetalScalar>(0.0);
    const uint32_t group_count = partial_count_for(active_count);
    for(uint32_t index = 0; index < group_count; ++index) {
      sum += partial[index];
    }
    return sum;
  }

  MetalScalar reduce_dot(id<MTLBuffer> left, id<MTLBuffer> right, const uint32_t count) {
    dispatch_1d(
        dot_partial_pipeline_,
        count,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:left offset:0 atIndex:0];
          [encoder setBuffer:right offset:0 atIndex:1];
          [encoder setBuffer:partial_buffer_ offset:0 atIndex:2];
          [encoder setBytes:&count length:sizeof(count) atIndex:3];
        });

    const MetalScalar* partial = buffer_data<MetalScalar>(partial_buffer_);
    MetalScalar sum = static_cast<MetalScalar>(0.0);
    const uint32_t group_count = partial_count_for(count);
    for(uint32_t index = 0; index < group_count; ++index) {
      sum += partial[index];
    }
    return sum;
  }

  double compute_max_abs_difference(id<MTLBuffer> left, id<MTLBuffer> right) {
    const uint32_t count = face_count_;
    dispatch_1d(
        max_abs_diff_pipeline_,
        count,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:left offset:0 atIndex:0];
          [encoder setBuffer:right offset:0 atIndex:1];
          [encoder setBuffer:partial_buffer_ offset:0 atIndex:2];
          [encoder setBytes:&count length:sizeof(count) atIndex:3];
        });

    const MetalScalar* partial = buffer_data<MetalScalar>(partial_buffer_);
    MetalScalar maximum = static_cast<MetalScalar>(0.0);
    const uint32_t group_count = partial_count_for(count);
    for(uint32_t index = 0; index < group_count; ++index) {
      maximum = std::max(maximum, partial[index]);
    }
    return static_cast<double>(maximum);
  }

  double compute_cfl(const std::array<id<MTLBuffer>, 3>& velocity) {
    dispatch_1d(
        cfl_partial_pipeline_,
        cell_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:velocity[0] offset:0 atIndex:0];
          [encoder setBuffer:velocity[1] offset:0 atIndex:1];
          [encoder setBuffer:velocity[2] offset:0 atIndex:2];
          [encoder setBuffer:partial_buffer_ offset:0 atIndex:3];
          [encoder setBytes:&params_ length:sizeof(params_) atIndex:4];
        });

    const MetalScalar* partial = buffer_data<MetalScalar>(partial_buffer_);
    MetalScalar maximum = static_cast<MetalScalar>(0.0);
    const uint32_t group_count = partial_count_for(cell_count_);
    for(uint32_t index = 0; index < group_count; ++index) {
      maximum = std::max(maximum, partial[index]);
    }
    return static_cast<double>(maximum);
  }

  double compute_divergence_l2(const std::array<id<MTLBuffer>, 3>& velocity) {
    dispatch_1d(
        divergence_partial_pipeline_,
        cell_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:velocity[0] offset:0 atIndex:0];
          [encoder setBuffer:velocity[1] offset:0 atIndex:1];
          [encoder setBuffer:velocity[2] offset:0 atIndex:2];
          [encoder setBuffer:partial_buffer_ offset:0 atIndex:3];
          [encoder setBytes:&params_ length:sizeof(params_) atIndex:4];
        });

    const MetalScalar* partial = buffer_data<MetalScalar>(partial_buffer_);
    MetalScalar sum = static_cast<MetalScalar>(0.0);
    const uint32_t group_count = partial_count_for(cell_count_);
    for(uint32_t index = 0; index < group_count; ++index) {
      sum += partial[index];
    }
    return std::sqrt(static_cast<double>(sum) / static_cast<double>(cell_count_));
  }

  void subtract_scalar_from_buffer(id<MTLBuffer> buffer,
                                   const MetalScalar value,
                                   const uint32_t count) {
    dispatch_1d(
        subtract_scalar_pipeline_,
        count,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:buffer offset:0 atIndex:0];
          [encoder setBytes:&value length:sizeof(value) atIndex:1];
          [encoder setBytes:&count length:sizeof(count) atIndex:2];
        });
  }

  void solve_face_system(id<MTLBuffer> solution) {
    copy_buffer(face_rhs_, solution, static_cast<NSUInteger>(face_count_) * sizeof(MetalScalar));

    dispatch_1d(
        factorized_operator_pipeline_,
        face_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:solution offset:0 atIndex:0];
          [encoder setBuffer:face_applied_ offset:0 atIndex:1];
          [encoder setBytes:&params_ length:sizeof(params_) atIndex:2];
        });

    dispatch_1d(
        init_residual_pipeline_,
        face_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:face_rhs_ offset:0 atIndex:0];
          [encoder setBuffer:face_applied_ offset:0 atIndex:1];
          [encoder setBuffer:face_residual_ offset:0 atIndex:2];
          [encoder setBuffer:face_direction_ offset:0 atIndex:3];
          [encoder setBytes:&face_count_ length:sizeof(face_count_) atIndex:4];
        });

    MetalScalar residual_dot = reduce_dot(face_residual_, face_residual_, face_count_);
    if(static_cast<double>(residual_dot) <= std::numeric_limits<double>::epsilon()) {
      return;
    }

    const MetalScalar initial_dot = residual_dot;
    for(int iteration = 0; iteration < 80; ++iteration) {
      dispatch_1d(
          factorized_operator_pipeline_,
          face_count_,
          [&](id<MTLComputeCommandEncoder> encoder) {
            [encoder setBuffer:face_direction_ offset:0 atIndex:0];
            [encoder setBuffer:face_applied_ offset:0 atIndex:1];
            [encoder setBytes:&params_ length:sizeof(params_) atIndex:2];
          });

      const MetalScalar search_dot = reduce_dot(face_direction_, face_applied_, face_count_);
      if(std::abs(static_cast<double>(search_dot)) <= std::numeric_limits<double>::epsilon()) {
        break;
      }

      const MetalScalar alpha = residual_dot / search_dot;
      dispatch_1d(
          update_solution_pipeline_,
          face_count_,
          [&](id<MTLComputeCommandEncoder> encoder) {
            [encoder setBuffer:solution offset:0 atIndex:0];
            [encoder setBuffer:face_residual_ offset:0 atIndex:1];
            [encoder setBuffer:face_direction_ offset:0 atIndex:2];
            [encoder setBuffer:face_applied_ offset:0 atIndex:3];
            [encoder setBytes:&alpha length:sizeof(alpha) atIndex:4];
            [encoder setBytes:&face_count_ length:sizeof(face_count_) atIndex:5];
          });

      const MetalScalar next_dot = reduce_dot(face_residual_, face_residual_, face_count_);
      const double relative = std::sqrt(static_cast<double>(next_dot / initial_dot));
      if(relative <= kCgTolerance) {
        break;
      }

      const MetalScalar beta = next_dot / residual_dot;
      dispatch_1d(
          update_direction_pipeline_,
          face_count_,
          [&](id<MTLComputeCommandEncoder> encoder) {
            [encoder setBuffer:face_direction_ offset:0 atIndex:0];
            [encoder setBuffer:face_residual_ offset:0 atIndex:1];
            [encoder setBytes:&beta length:sizeof(beta) atIndex:2];
            [encoder setBytes:&face_count_ length:sizeof(face_count_) atIndex:3];
          });
      residual_dot = next_dot;
    }
  }

  void build_pressure_rhs() {
    dispatch_1d(
        pressure_rhs_pipeline_,
        cell_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:next_velocity_[0] offset:0 atIndex:0];
          [encoder setBuffer:next_velocity_[1] offset:0 atIndex:1];
          [encoder setBuffer:next_velocity_[2] offset:0 atIndex:2];
          [encoder setBuffer:pressure_rhs_ offset:0 atIndex:3];
          [encoder setBytes:&params_ length:sizeof(params_) atIndex:4];
        });
  }

  PressureSolveStats solve_pressure_system() {
    zero_buffer(pressure_correction_);
    copy_buffer(pressure_rhs_,
                pressure_residual_,
                static_cast<NSUInteger>(cell_count_) * sizeof(MetalScalar));
    copy_buffer(pressure_rhs_,
                pressure_direction_,
                static_cast<NSUInteger>(cell_count_) * sizeof(MetalScalar));

    MetalScalar residual_dot = reduce_dot(pressure_residual_, pressure_residual_, cell_count_);
    if(static_cast<double>(residual_dot) <= std::numeric_limits<double>::epsilon()) {
      return PressureSolveStats{};
    }

    const MetalScalar initial_dot = residual_dot;
    PressureSolveStats stats{};
    stats.relative_residual = 1.0;

    const double target_relative = std::max(config_.poisson_tolerance, kCgTolerance);
    for(int iteration = 0; iteration < config_.poisson_max_iterations; ++iteration) {
      dispatch_1d(
          poisson_operator_pipeline_,
          cell_count_,
          [&](id<MTLComputeCommandEncoder> encoder) {
            [encoder setBuffer:pressure_direction_ offset:0 atIndex:0];
            [encoder setBuffer:pressure_applied_ offset:0 atIndex:1];
            [encoder setBytes:&params_ length:sizeof(params_) atIndex:2];
          });

      const MetalScalar search_dot =
          reduce_dot(pressure_direction_, pressure_applied_, cell_count_);
      if(std::abs(static_cast<double>(search_dot)) <= std::numeric_limits<double>::epsilon()) {
        break;
      }

      const MetalScalar alpha = residual_dot / search_dot;
      dispatch_1d(
          update_solution_pipeline_,
          cell_count_,
          [&](id<MTLComputeCommandEncoder> encoder) {
            [encoder setBuffer:pressure_correction_ offset:0 atIndex:0];
            [encoder setBuffer:pressure_residual_ offset:0 atIndex:1];
            [encoder setBuffer:pressure_direction_ offset:0 atIndex:2];
            [encoder setBuffer:pressure_applied_ offset:0 atIndex:3];
            [encoder setBytes:&alpha length:sizeof(alpha) atIndex:4];
            [encoder setBytes:&cell_count_ length:sizeof(cell_count_) atIndex:5];
          });

      const MetalScalar next_dot = reduce_dot(pressure_residual_, pressure_residual_, cell_count_);
      stats.iterations = iteration + 1;
      stats.relative_residual = std::sqrt(static_cast<double>(next_dot / initial_dot));
      if(stats.relative_residual <= target_relative) {
        break;
      }

      const MetalScalar beta = next_dot / residual_dot;
      dispatch_1d(
          update_direction_pipeline_,
          cell_count_,
          [&](id<MTLComputeCommandEncoder> encoder) {
            [encoder setBuffer:pressure_direction_ offset:0 atIndex:0];
            [encoder setBuffer:pressure_residual_ offset:0 atIndex:1];
            [encoder setBytes:&beta length:sizeof(beta) atIndex:2];
            [encoder setBytes:&cell_count_ length:sizeof(cell_count_) atIndex:3];
          });
      residual_dot = next_dot;
    }

    const MetalScalar mean = reduce_sum(pressure_correction_) / static_cast<MetalScalar>(cell_count_);
    if(std::abs(static_cast<double>(mean)) > 0.0) {
      subtract_scalar_from_buffer(pressure_correction_, mean, cell_count_);
    }
    return stats;
  }

  void correct_velocity(const uint32_t component) {
    id<MTLComputePipelineState> pipeline = correct_u_pipeline_;
    switch(component) {
      case 0u:
        pipeline = correct_u_pipeline_;
        break;
      case 1u:
        pipeline = correct_v_pipeline_;
        break;
      case 2u:
        pipeline = correct_w_pipeline_;
        break;
    }

    dispatch_1d(
        pipeline,
        face_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:next_velocity_[component] offset:0 atIndex:0];
          [encoder setBuffer:pressure_correction_ offset:0 atIndex:1];
          [encoder setBuffer:next_velocity_[component] offset:0 atIndex:2];
          [encoder setBytes:&params_ length:sizeof(params_) atIndex:3];
        });
  }

  void add_pressure_correction() {
    dispatch_1d(
        add_pressure_pipeline_,
        cell_count_,
        [&](id<MTLComputeCommandEncoder> encoder) {
          [encoder setBuffer:pressure_total_ offset:0 atIndex:0];
          [encoder setBuffer:pressure_correction_ offset:0 atIndex:1];
          [encoder setBytes:&cell_count_ length:sizeof(cell_count_) atIndex:2];
        });
  }

  void unpack_periodic_face_buffer(id<MTLBuffer> buffer, FaceField& field) {
    const MetalScalar* source = buffer_data<MetalScalar>(buffer);
    const IndexRange3D active = field.layout().active_range();
    const int axis = axis_index(field.normal_axis());

    for(int k = 0; k < grid_.nz; ++k) {
      for(int j = 0; j < grid_.ny; ++j) {
        for(int i = 0; i < grid_.nx; ++i) {
          const std::size_t index = static_cast<std::size_t>(i) +
                                    static_cast<std::size_t>(grid_.nx) *
                                        (static_cast<std::size_t>(j) +
                                         static_cast<std::size_t>(grid_.ny) *
                                             static_cast<std::size_t>(k));
          const int si = active.i_begin + i;
          const int sj = active.j_begin + j;
          const int sk = active.k_begin + k;
          field(si, sj, sk) = static_cast<double>(source[index]);
        }
      }
    }

    if(axis == axis_index(Axis::x)) {
      for(int k = active.k_begin; k < active.k_end; ++k) {
        for(int j = active.j_begin; j < active.j_end; ++j) {
          field(active.i_end - 1, j, k) = field(active.i_begin, j, k);
        }
      }
    } else if(axis == axis_index(Axis::y)) {
      for(int k = active.k_begin; k < active.k_end; ++k) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          field(i, active.j_end - 1, k) = field(i, active.j_begin, k);
        }
      }
    } else {
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          field(i, j, active.k_end - 1) = field(i, j, active.k_begin);
        }
      }
    }
  }

  void unpack_pressure_buffer(id<MTLBuffer> buffer, PressureField& field) {
    const MetalScalar* source = buffer_data<MetalScalar>(buffer);
    const IndexRange3D active = field.layout().active_range();
    for(int k = 0; k < grid_.nz; ++k) {
      for(int j = 0; j < grid_.ny; ++j) {
        for(int i = 0; i < grid_.nx; ++i) {
          const std::size_t index = static_cast<std::size_t>(i) +
                                    static_cast<std::size_t>(grid_.nx) *
                                        (static_cast<std::size_t>(j) +
                                         static_cast<std::size_t>(grid_.ny) *
                                             static_cast<std::size_t>(k));
          field(active.i_begin + i, active.j_begin + j, active.k_begin + k) =
              static_cast<double>(source[index]);
        }
      }
    }
  }

  TaylorGreenState unpack_state(const TaylorGreenStepMetrics& metrics,
                                const bool has_previous_advection) {
    TaylorGreenState state{grid_};
    unpack_periodic_face_buffer(current_velocity_[0], state.velocity.x);
    unpack_periodic_face_buffer(current_velocity_[1], state.velocity.y);
    unpack_periodic_face_buffer(current_velocity_[2], state.velocity.z);
    unpack_pressure_buffer(pressure_total_, state.pressure_total);
    state.metrics = metrics;
    state.has_previous_advection = has_previous_advection;
    return state;
  }

  TaylorGreenConfig config_;
  Grid grid_;
  int step_count_ = 0;
  double dt_ = 0.0;
  SolverParams params_{};
  uint32_t face_count_ = 0;
  uint32_t cell_count_ = 0;
  uint32_t partial_count_ = 0;

  id<MTLDevice> device_ = nil;
  id<MTLCommandQueue> command_queue_ = nil;
  id<MTLLibrary> library_ = nil;

  id<MTLComputePipelineState> init_pipeline_ = nil;
  id<MTLComputePipelineState> assemble_u_pipeline_ = nil;
  id<MTLComputePipelineState> assemble_v_pipeline_ = nil;
  id<MTLComputePipelineState> assemble_w_pipeline_ = nil;
  id<MTLComputePipelineState> factorized_operator_pipeline_ = nil;
  id<MTLComputePipelineState> init_residual_pipeline_ = nil;
  id<MTLComputePipelineState> update_solution_pipeline_ = nil;
  id<MTLComputePipelineState> update_direction_pipeline_ = nil;
  id<MTLComputePipelineState> subtract_scalar_pipeline_ = nil;
  id<MTLComputePipelineState> pressure_rhs_pipeline_ = nil;
  id<MTLComputePipelineState> poisson_operator_pipeline_ = nil;
  id<MTLComputePipelineState> correct_u_pipeline_ = nil;
  id<MTLComputePipelineState> correct_v_pipeline_ = nil;
  id<MTLComputePipelineState> correct_w_pipeline_ = nil;
  id<MTLComputePipelineState> add_pressure_pipeline_ = nil;
  id<MTLComputePipelineState> dot_partial_pipeline_ = nil;
  id<MTLComputePipelineState> sum_partial_pipeline_ = nil;
  id<MTLComputePipelineState> max_abs_diff_pipeline_ = nil;
  id<MTLComputePipelineState> cfl_partial_pipeline_ = nil;
  id<MTLComputePipelineState> divergence_partial_pipeline_ = nil;

  std::array<id<MTLBuffer>, 3> current_velocity_{};
  std::array<id<MTLBuffer>, 3> next_velocity_{};
  std::array<id<MTLBuffer>, 3> previous_advection_{};
  std::array<id<MTLBuffer>, 3> current_advection_{};

  id<MTLBuffer> pressure_total_ = nil;
  id<MTLBuffer> pressure_rhs_ = nil;
  id<MTLBuffer> pressure_correction_ = nil;
  id<MTLBuffer> pressure_residual_ = nil;
  id<MTLBuffer> pressure_direction_ = nil;
  id<MTLBuffer> pressure_applied_ = nil;

  id<MTLBuffer> face_rhs_ = nil;
  id<MTLBuffer> face_residual_ = nil;
  id<MTLBuffer> face_direction_ = nil;
  id<MTLBuffer> face_applied_ = nil;

  id<MTLBuffer> partial_buffer_ = nil;
};

}  // namespace

TaylorGreenMetalRun run_taylor_green(const TaylorGreenConfig& config) {
  MetalRunner runner{config};
  return runner.run();
}

}  // namespace solver::metal
