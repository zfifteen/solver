#!/bin/zsh

set -euo pipefail

dry_run=0
if [[ "${1:-}" == "--dry-run" ]]; then
  dry_run=1
  shift
fi

build_dir="${1:-build/benchmark}"
output_dir="${2:-profiling/latest}"

script_dir="$(cd "$(dirname "$0")" && pwd)"
root_dir="$(cd "${script_dir}/.." && pwd)"

mkdir -p "${root_dir}/${output_dir}"

timestamp="$(date +"%Y%m%d-%H%M%S")"
log_path="${root_dir}/${output_dir}/run_profile_suite-${timestamp}.log"
pid_path="${root_dir}/${output_dir}/run_profile_suite.pid"

command=(
  "./profiling/run_profile_suite.py"
  "--build-dir"
  "${build_dir}"
  "--output-dir"
  "${output_dir}"
)

if (( dry_run )); then
  echo "cwd: ${root_dir}"
  echo "log: ${log_path}"
  echo "pid_file: ${pid_path}"
  printf "command:"
  printf " %q" "${command[@]}"
  printf "\n"
  exit 0
fi

cd "${root_dir}"
nohup "${command[@]}" >"${log_path}" 2>&1 < /dev/null &
pid=$!

echo "${pid}" > "${pid_path}"

echo "started profile suite"
echo "pid: ${pid}"
echo "log: ${log_path}"
echo "pid_file: ${pid_path}"
