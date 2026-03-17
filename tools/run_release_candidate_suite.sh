#!/bin/zsh

set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
root_dir="$(cd "${script_dir}/.." && pwd)"

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/solver-rc-suite.XXXXXX")"
checkpoint_path="${tmp_dir}/solver.chk"

cleanup() {
  rm -rf "${tmp_dir}"
}
trap cleanup EXIT

cd "${root_dir}"

echo "[rc] configure deterministic"
cmake --preset deterministic

echo "[rc] build deterministic"
cmake --build --preset deterministic -j 3

echo "[rc] run deterministic tests"
ctest --preset deterministic

echo "[rc] run deterministic validation harness"
./validation/run_validation_suite.py --build-dir build/deterministic --output-dir validation/latest

echo "[rc] write cavity checkpoint smoke"
build/deterministic/tools/solver_cavity benchmarks/lid_driven_cavity_smoke.cfg \
  --steps 6 \
  --checkpoint-out "${checkpoint_path}"

echo "[rc] restart cavity checkpoint smoke"
build/deterministic/tools/solver_cavity benchmarks/lid_driven_cavity_smoke.cfg \
  --restart "${checkpoint_path}" \
  --steps 6

echo "[rc] release-candidate deterministic suite completed"
