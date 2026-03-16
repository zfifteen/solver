#include "io/checkpoint.hpp"

#include "core/runtime.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace solver::io {

namespace {

using ByteBuffer = std::vector<std::uint8_t>;

constexpr std::array<char, 8> kMagic{'S', 'L', 'V', 'C', 'H', 'K', 'P', '1'};
constexpr std::uint32_t kFormatVersion = 1;

bool same_grid(const Grid& left, const Grid& right) {
  return left.nx == right.nx && left.ny == right.ny && left.nz == right.nz &&
         left.dx == right.dx && left.dy == right.dy && left.dz == right.dz &&
         left.ghost_layers == right.ghost_layers;
}

std::uint64_t fnv1a_64(const std::uint8_t* data, const std::size_t size) {
  std::uint64_t hash = 14695981039346656037ull;
  for(std::size_t index = 0; index < size; ++index) {
    hash ^= static_cast<std::uint64_t>(data[index]);
    hash *= 1099511628211ull;
  }
  return hash;
}

std::uint64_t hash_string(const std::string& text) {
  return fnv1a_64(reinterpret_cast<const std::uint8_t*>(text.data()), text.size());
}

std::uint64_t build_hash() {
  const BuildInfo build_info = get_build_info();
  const std::string description = build_info.project_name + "|" + build_info.project_version + "|" +
                                  build_info.compiler_id + "|" + build_info.compiler_version + "|" +
                                  build_info.build_profile + "|" +
                                  (build_info.fast_math_enabled ? "fast-math" : "strict");
  return hash_string(description);
}

std::uint64_t configuration_hash(const LidDrivenCavityConfig& config) {
  return hash_string(describe(config));
}

void append_u32(ByteBuffer& buffer, const std::uint32_t value) {
  for(int shift = 0; shift < 32; shift += 8) {
    buffer.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffu));
  }
}

void append_u64(ByteBuffer& buffer, const std::uint64_t value) {
  for(int shift = 0; shift < 64; shift += 8) {
    buffer.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffu));
  }
}

void append_double(ByteBuffer& buffer, const double value) {
  append_u64(buffer, std::bit_cast<std::uint64_t>(value));
}

void append_tag(ByteBuffer& buffer, const std::array<char, 4>& tag) {
  buffer.insert(buffer.end(), tag.begin(), tag.end());
}

std::uint32_t read_u32(const ByteBuffer& buffer, std::size_t& cursor) {
  if(cursor + 4 > buffer.size()) {
    throw std::runtime_error("checkpoint truncated while reading u32");
  }

  std::uint32_t value = 0;
  for(int shift = 0; shift < 32; shift += 8) {
    value |= static_cast<std::uint32_t>(buffer[cursor++]) << shift;
  }
  return value;
}

std::uint64_t read_u64(const ByteBuffer& buffer, std::size_t& cursor) {
  if(cursor + 8 > buffer.size()) {
    throw std::runtime_error("checkpoint truncated while reading u64");
  }

  std::uint64_t value = 0;
  for(int shift = 0; shift < 64; shift += 8) {
    value |= static_cast<std::uint64_t>(buffer[cursor++]) << shift;
  }
  return value;
}

double read_double(const ByteBuffer& buffer, std::size_t& cursor) {
  return std::bit_cast<double>(read_u64(buffer, cursor));
}

ByteBuffer serialize_metadata(const LidDrivenCavityConfig& config,
                              const LidDrivenCavityState& state) {
  ByteBuffer bytes;
  append_u32(bytes, static_cast<std::uint32_t>(sizeof(double)));
  append_u64(bytes, build_hash());
  append_u64(bytes, configuration_hash(config));
  append_u32(bytes, static_cast<std::uint32_t>(state.has_previous_advection ? 1 : 0));
  append_u32(bytes, static_cast<std::uint32_t>(state.metrics.step));
  append_double(bytes, state.metrics.time);
  append_double(bytes, state.metrics.dt);
  append_double(bytes, state.metrics.max_cfl);
  append_double(bytes, state.metrics.max_velocity_change);
  append_double(bytes, state.metrics.divergence_l2);
  append_u32(bytes, static_cast<std::uint32_t>(state.metrics.pressure_iterations));
  append_double(bytes, state.metrics.pressure_relative_residual);
  return bytes;
}

ByteBuffer serialize_config(const LidDrivenCavityConfig& config) {
  ByteBuffer bytes;
  append_u32(bytes, static_cast<std::uint32_t>(config.nx));
  append_u32(bytes, static_cast<std::uint32_t>(config.ny));
  append_double(bytes, config.reynolds);
  append_double(bytes, config.lid_velocity);
  append_double(bytes, config.cfl_limit);
  append_u32(bytes, static_cast<std::uint32_t>(config.max_steps));
  append_u32(bytes, static_cast<std::uint32_t>(config.min_steps));
  append_double(bytes, config.steady_tolerance);
  append_u32(bytes, static_cast<std::uint32_t>(config.poisson_max_iterations));
  append_double(bytes, config.poisson_tolerance);
  append_u32(bytes, static_cast<std::uint32_t>(config.validate_reference ? 1 : 0));
  append_u32(bytes, static_cast<std::uint32_t>(config.advection.scheme));
  append_u32(bytes, static_cast<std::uint32_t>(config.advection.limiter));
  return bytes;
}

ByteBuffer serialize_grid(const Grid& grid) {
  ByteBuffer bytes;
  append_u32(bytes, static_cast<std::uint32_t>(grid.nx));
  append_u32(bytes, static_cast<std::uint32_t>(grid.ny));
  append_u32(bytes, static_cast<std::uint32_t>(grid.nz));
  append_double(bytes, grid.dx);
  append_double(bytes, grid.dy);
  append_double(bytes, grid.dz);
  append_u32(bytes, static_cast<std::uint32_t>(grid.ghost_layers));
  return bytes;
}

ByteBuffer serialize_field(const StructuredField& field) {
  ByteBuffer bytes;
  const Extent3D storage = field.layout().storage_extent();
  append_u32(bytes, static_cast<std::uint32_t>(storage.nx));
  append_u32(bytes, static_cast<std::uint32_t>(storage.ny));
  append_u32(bytes, static_cast<std::uint32_t>(storage.nz));

  for(int k = 0; k < storage.nz; ++k) {
    for(int j = 0; j < storage.ny; ++j) {
      for(int i = 0; i < storage.nx; ++i) {
        append_double(bytes, field(i, j, k));
      }
    }
  }

  return bytes;
}

void append_section(ByteBuffer& payload,
                    const std::array<char, 4>& tag,
                    const ByteBuffer& section_bytes) {
  append_tag(payload, tag);
  append_u64(payload, static_cast<std::uint64_t>(section_bytes.size()));
  payload.insert(payload.end(), section_bytes.begin(), section_bytes.end());
}

void validate_serialized_config(const ByteBuffer& bytes,
                                const LidDrivenCavityConfig& expected_config) {
  std::size_t cursor = 0;
  LidDrivenCavityConfig loaded{};
  loaded.nx = static_cast<int>(read_u32(bytes, cursor));
  loaded.ny = static_cast<int>(read_u32(bytes, cursor));
  loaded.reynolds = read_double(bytes, cursor);
  loaded.lid_velocity = read_double(bytes, cursor);
  loaded.cfl_limit = read_double(bytes, cursor);
  loaded.max_steps = static_cast<int>(read_u32(bytes, cursor));
  loaded.min_steps = static_cast<int>(read_u32(bytes, cursor));
  loaded.steady_tolerance = read_double(bytes, cursor);
  loaded.poisson_max_iterations = static_cast<int>(read_u32(bytes, cursor));
  loaded.poisson_tolerance = read_double(bytes, cursor);
  loaded.validate_reference = read_u32(bytes, cursor) != 0u;
  loaded.advection.scheme = static_cast<AdvectionScheme>(read_u32(bytes, cursor));
  loaded.advection.limiter = static_cast<FluxLimiter>(read_u32(bytes, cursor));

  if(cursor != bytes.size()) {
    throw std::runtime_error("checkpoint config section contains trailing bytes");
  }

  if(loaded.nx != expected_config.nx || loaded.ny != expected_config.ny ||
     loaded.reynolds != expected_config.reynolds ||
     loaded.lid_velocity != expected_config.lid_velocity ||
     loaded.cfl_limit != expected_config.cfl_limit ||
     loaded.max_steps != expected_config.max_steps ||
     loaded.min_steps != expected_config.min_steps ||
     loaded.steady_tolerance != expected_config.steady_tolerance ||
     loaded.poisson_max_iterations != expected_config.poisson_max_iterations ||
     loaded.poisson_tolerance != expected_config.poisson_tolerance ||
     loaded.validate_reference != expected_config.validate_reference ||
     loaded.advection.scheme != expected_config.advection.scheme ||
     loaded.advection.limiter != expected_config.advection.limiter) {
    throw std::runtime_error("checkpoint config does not match the requested run configuration");
  }
}

void validate_serialized_grid(const ByteBuffer& bytes, const Grid& expected_grid) {
  std::size_t cursor = 0;
  const Grid loaded{
      static_cast<int>(read_u32(bytes, cursor)),
      static_cast<int>(read_u32(bytes, cursor)),
      static_cast<int>(read_u32(bytes, cursor)),
      read_double(bytes, cursor),
      read_double(bytes, cursor),
      read_double(bytes, cursor),
      static_cast<int>(read_u32(bytes, cursor)),
  };

  if(cursor != bytes.size()) {
    throw std::runtime_error("checkpoint grid section contains trailing bytes");
  }

  if(!same_grid(loaded, expected_grid)) {
    throw std::runtime_error("checkpoint grid does not match the requested run configuration");
  }
}

void deserialize_metadata(const ByteBuffer& bytes,
                          const LidDrivenCavityConfig& config,
                          LidDrivenCavityCheckpoint& checkpoint) {
  std::size_t cursor = 0;
  checkpoint.metadata.scalar_bytes = read_u32(bytes, cursor);
  checkpoint.metadata.build_hash = read_u64(bytes, cursor);
  checkpoint.metadata.configuration_hash = read_u64(bytes, cursor);
  checkpoint.state.has_previous_advection = read_u32(bytes, cursor) != 0u;
  checkpoint.state.metrics.step = static_cast<int>(read_u32(bytes, cursor));
  checkpoint.state.metrics.time = read_double(bytes, cursor);
  checkpoint.state.metrics.dt = read_double(bytes, cursor);
  checkpoint.state.metrics.max_cfl = read_double(bytes, cursor);
  checkpoint.state.metrics.max_velocity_change = read_double(bytes, cursor);
  checkpoint.state.metrics.divergence_l2 = read_double(bytes, cursor);
  checkpoint.state.metrics.pressure_iterations = static_cast<int>(read_u32(bytes, cursor));
  checkpoint.state.metrics.pressure_relative_residual = read_double(bytes, cursor);

  if(cursor != bytes.size()) {
    throw std::runtime_error("checkpoint metadata section contains trailing bytes");
  }

  if(checkpoint.metadata.scalar_bytes != sizeof(double)) {
    throw std::runtime_error("checkpoint precision does not match the current solver precision");
  }
  if(checkpoint.metadata.build_hash != build_hash()) {
    throw std::runtime_error("checkpoint build hash does not match the current executable");
  }
  if(checkpoint.metadata.configuration_hash != configuration_hash(config)) {
    throw std::runtime_error("checkpoint configuration hash does not match the requested config");
  }
}

void deserialize_field(const ByteBuffer& bytes, StructuredField& field) {
  std::size_t cursor = 0;
  const Extent3D expected = field.layout().storage_extent();
  const int nx = static_cast<int>(read_u32(bytes, cursor));
  const int ny = static_cast<int>(read_u32(bytes, cursor));
  const int nz = static_cast<int>(read_u32(bytes, cursor));
  if(nx != expected.nx || ny != expected.ny || nz != expected.nz) {
    throw std::runtime_error("checkpoint field extent does not match the current layout");
  }

  for(int k = 0; k < expected.nz; ++k) {
    for(int j = 0; j < expected.ny; ++j) {
      for(int i = 0; i < expected.nx; ++i) {
        field(i, j, k) = read_double(bytes, cursor);
      }
    }
  }

  if(cursor != bytes.size()) {
    throw std::runtime_error("checkpoint field section contains trailing bytes");
  }
}

ByteBuffer read_file_bytes(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if(!input.is_open()) {
    throw std::runtime_error("unable to open checkpoint: " + path);
  }

  input.seekg(0, std::ios::end);
  const std::streamsize size = input.tellg();
  input.seekg(0, std::ios::beg);
  if(size < 0) {
    throw std::runtime_error("unable to read checkpoint size: " + path);
  }

  ByteBuffer bytes(static_cast<std::size_t>(size), 0u);
  if(!input.read(reinterpret_cast<char*>(bytes.data()), size)) {
    throw std::runtime_error("unable to read checkpoint payload: " + path);
  }
  return bytes;
}

}  // namespace

void write_lid_driven_cavity_checkpoint(const std::string& path,
                                        const LidDrivenCavityConfig& config,
                                        const LidDrivenCavityState& state) {
  const Grid expected_grid = make_lid_driven_cavity_grid(config);
  if(!same_grid(state.grid, expected_grid)) {
    throw std::invalid_argument("checkpoint state grid does not match the cavity config");
  }

  ByteBuffer payload;
  append_section(payload, {'M', 'E', 'T', 'A'}, serialize_metadata(config, state));
  append_section(payload, {'C', 'O', 'N', 'F'}, serialize_config(config));
  append_section(payload, {'G', 'R', 'I', 'D'}, serialize_grid(state.grid));
  append_section(payload, {'V', 'E', 'L', 'X'}, serialize_field(state.velocity.x));
  append_section(payload, {'V', 'E', 'L', 'Y'}, serialize_field(state.velocity.y));
  append_section(payload, {'V', 'E', 'L', 'Z'}, serialize_field(state.velocity.z));
  append_section(payload, {'A', 'D', 'V', 'X'}, serialize_field(state.advection_previous.x));
  append_section(payload, {'A', 'D', 'V', 'Y'}, serialize_field(state.advection_previous.y));
  append_section(payload, {'A', 'D', 'V', 'Z'}, serialize_field(state.advection_previous.z));
  append_section(payload, {'P', 'R', 'S', 'T'}, serialize_field(state.pressure_total));

  const std::uint64_t checksum = fnv1a_64(payload.data(), payload.size());
  std::ofstream output(path, std::ios::binary);
  if(!output.is_open()) {
    throw std::runtime_error("unable to open checkpoint for writing: " + path);
  }

  output.write(kMagic.data(), static_cast<std::streamsize>(kMagic.size()));
  ByteBuffer header;
  append_u32(header, kFormatVersion);
  append_u64(header, static_cast<std::uint64_t>(payload.size()));
  append_u64(header, checksum);
  output.write(reinterpret_cast<const char*>(header.data()),
               static_cast<std::streamsize>(header.size()));
  output.write(reinterpret_cast<const char*>(payload.data()),
               static_cast<std::streamsize>(payload.size()));
  if(!output.good()) {
    throw std::runtime_error("failed while writing checkpoint: " + path);
  }
}

LidDrivenCavityCheckpoint load_lid_driven_cavity_checkpoint(const std::string& path,
                                                            const LidDrivenCavityConfig& expected_config) {
  const Grid expected_grid = make_lid_driven_cavity_grid(expected_config);
  ByteBuffer bytes = read_file_bytes(path);
  const std::size_t minimum_size = kMagic.size() + 4u + 8u + 8u;
  if(bytes.size() < minimum_size) {
    throw std::runtime_error("checkpoint file is too small: " + path);
  }

  std::size_t cursor = 0;
  for(const char expected_byte : kMagic) {
    if(static_cast<char>(bytes[cursor++]) != expected_byte) {
      throw std::runtime_error("checkpoint magic mismatch: " + path);
    }
  }

  const std::uint32_t format_version = read_u32(bytes, cursor);
  if(format_version != kFormatVersion) {
    throw std::runtime_error("unsupported checkpoint format version: " + std::to_string(format_version));
  }

  const std::uint64_t payload_size = read_u64(bytes, cursor);
  const std::uint64_t checksum = read_u64(bytes, cursor);
  if(cursor + payload_size != bytes.size()) {
    throw std::runtime_error("checkpoint payload length does not match the file size");
  }

  const std::uint64_t computed_checksum =
      fnv1a_64(bytes.data() + cursor, static_cast<std::size_t>(payload_size));
  if(computed_checksum != checksum) {
    throw std::runtime_error("checkpoint checksum verification failed");
  }

  LidDrivenCavityCheckpoint checkpoint{expected_grid};
  checkpoint.metadata.format_version = format_version;
  checkpoint.metadata.checksum = checksum;

  bool has_metadata = false;
  bool has_config = false;
  bool has_grid = false;
  bool has_velx = false;
  bool has_vely = false;
  bool has_velz = false;
  bool has_advx = false;
  bool has_advy = false;
  bool has_advz = false;
  bool has_pressure = false;

  while(cursor < bytes.size()) {
    if(cursor + 4u + 8u > bytes.size()) {
      throw std::runtime_error("checkpoint section header is truncated");
    }

    const std::array<char, 4> tag{
        static_cast<char>(bytes[cursor]),
        static_cast<char>(bytes[cursor + 1]),
        static_cast<char>(bytes[cursor + 2]),
        static_cast<char>(bytes[cursor + 3]),
    };
    cursor += 4u;
    const std::uint64_t section_size = read_u64(bytes, cursor);
    if(cursor + section_size > bytes.size()) {
      throw std::runtime_error("checkpoint section payload is truncated");
    }

    const ByteBuffer section(bytes.begin() + static_cast<std::ptrdiff_t>(cursor),
                             bytes.begin() + static_cast<std::ptrdiff_t>(cursor + section_size));
    cursor += static_cast<std::size_t>(section_size);

    if(tag == std::array<char, 4>{'M', 'E', 'T', 'A'}) {
      deserialize_metadata(section, expected_config, checkpoint);
      has_metadata = true;
    } else if(tag == std::array<char, 4>{'C', 'O', 'N', 'F'}) {
      validate_serialized_config(section, expected_config);
      has_config = true;
    } else if(tag == std::array<char, 4>{'G', 'R', 'I', 'D'}) {
      validate_serialized_grid(section, expected_grid);
      has_grid = true;
    } else if(tag == std::array<char, 4>{'V', 'E', 'L', 'X'}) {
      deserialize_field(section, checkpoint.state.velocity.x);
      has_velx = true;
    } else if(tag == std::array<char, 4>{'V', 'E', 'L', 'Y'}) {
      deserialize_field(section, checkpoint.state.velocity.y);
      has_vely = true;
    } else if(tag == std::array<char, 4>{'V', 'E', 'L', 'Z'}) {
      deserialize_field(section, checkpoint.state.velocity.z);
      has_velz = true;
    } else if(tag == std::array<char, 4>{'A', 'D', 'V', 'X'}) {
      deserialize_field(section, checkpoint.state.advection_previous.x);
      has_advx = true;
    } else if(tag == std::array<char, 4>{'A', 'D', 'V', 'Y'}) {
      deserialize_field(section, checkpoint.state.advection_previous.y);
      has_advy = true;
    } else if(tag == std::array<char, 4>{'A', 'D', 'V', 'Z'}) {
      deserialize_field(section, checkpoint.state.advection_previous.z);
      has_advz = true;
    } else if(tag == std::array<char, 4>{'P', 'R', 'S', 'T'}) {
      deserialize_field(section, checkpoint.state.pressure_total);
      has_pressure = true;
    }
  }

  if(!has_metadata || !has_config || !has_grid || !has_velx || !has_vely || !has_velz ||
     !has_advx || !has_advy || !has_advz || !has_pressure) {
    throw std::runtime_error("checkpoint is missing one or more required sections");
  }

  return checkpoint;
}

}  // namespace solver::io
