#!/bin/zsh

set -euo pipefail

build_dir="${1:-build/deterministic}"
binary_path="${build_dir}/tools/solver_example"

if [[ ! -x "${binary_path}" ]]; then
  echo "expected executable at ${binary_path}" >&2
  echo "build it first with: cmake --build ${build_dir}" >&2
  exit 1
fi

echo "profiling ${binary_path}"
/usr/bin/time -lp "${binary_path}"

