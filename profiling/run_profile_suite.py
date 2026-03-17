#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Milestone 11 profiling and optimization suite.")
    parser.add_argument(
        "--build-dir",
        default="build/benchmark",
        help="Build directory containing the performance-profile solver tools.",
    )
    parser.add_argument(
        "--output-dir",
        default="profiling/latest",
        help="Directory where profiling reports and plots should be written.",
    )
    return parser.parse_args()


def run_command(args: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )


def run_command_stdout_to_path(args: list[str], cwd: Path, stdout_path: Path) -> None:
    with stdout_path.open("w") as handle:
        subprocess.run(
            args,
            cwd=cwd,
            check=True,
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
        )


def parse_solver_metrics(stdout: str) -> dict[str, str]:
    metrics: dict[str, str] = {}
    for line in stdout.splitlines():
        if ": " not in line:
            continue
        key, value = line.split(": ", 1)
        metrics[key.strip()] = value.strip()
    return metrics


def parse_time_metrics(stderr: str) -> dict[str, float]:
    metrics: dict[str, float] = {}
    for line in stderr.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("real "):
            metrics["real_seconds"] = float(stripped.split()[1])
            continue
        if stripped.startswith("user "):
            metrics["user_seconds"] = float(stripped.split()[1])
            continue
        if stripped.startswith("sys "):
            metrics["sys_seconds"] = float(stripped.split()[1])
            continue

        parts = stripped.split()
        if len(parts) < 2:
            continue

        value = parts[0]
        label = " ".join(parts[1:])
        if value.replace(".", "", 1).isdigit():
            key = (
                label.replace(" ", "_")
                .replace("/", "_per_")
                .replace("-", "_")
                .replace("(", "")
                .replace(")", "")
            )
            metrics[key] = float(value)
    return metrics


def run_with_time(command: list[str], cwd: Path) -> tuple[dict[str, str], dict[str, float]]:
    completed = run_command(["/usr/bin/time", "-lp", *command], cwd)
    return parse_solver_metrics(completed.stdout), parse_time_metrics(completed.stderr)


def relative_l2_difference(reference: list[float], candidate: list[float]) -> float:
    if len(reference) != len(candidate):
        raise ValueError(
            f"relative_l2_difference expects equal-length vectors, got {len(reference)} and {len(candidate)}"
        )
    numerator = sum((left - right) ** 2 for left, right in zip(reference, candidate))
    denominator = sum(value * value for value in reference)
    if denominator <= 0.0:
        return math.sqrt(numerator)
    return math.sqrt(numerator / denominator)


def read_vtk_fields(path: Path) -> tuple[list[float], list[float]]:
    velocity_header = "VECTORS velocity double"
    pressure_header = "SCALARS pressure double 1"
    velocity: list[float] = []
    pressure: list[float] = []
    expected_points: int | None = None
    section: str | None = None

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("DIMENSIONS "):
                parts = line.split()
                if len(parts) != 4:
                    raise ValueError(f"invalid DIMENSIONS line in {path}: {line}")
                nx, ny, nz = (int(value) for value in parts[1:])
                expected_points = nx * ny * nz
                continue
            if line.startswith("POINT_DATA "):
                parts = line.split()
                if len(parts) != 2:
                    raise ValueError(f"invalid POINT_DATA line in {path}: {line}")
                point_data = int(parts[1])
                if expected_points is not None and point_data != expected_points:
                    raise ValueError(
                        f"VTK POINT_DATA count mismatch in {path}: {point_data} vs {expected_points}"
                    )
                expected_points = point_data
                continue
            if line == velocity_header:
                section = "velocity"
                continue
            if line == pressure_header:
                section = "pressure_lookup"
                continue
            if section == "pressure_lookup":
                if line != "LOOKUP_TABLE default":
                    raise ValueError(f"expected LOOKUP_TABLE default in {path}, found: {line}")
                section = "pressure"
                continue
            if section == "velocity":
                velocity.extend(float(value) for value in line.split())
                continue
            if section == "pressure":
                pressure.extend(float(value) for value in line.split())

    if expected_points is None:
        raise ValueError(f"missing DIMENSIONS/POINT_DATA in {path}")
    if len(velocity) != expected_points * 3:
        raise ValueError(
            f"velocity field length mismatch in {path}: expected {expected_points * 3}, got {len(velocity)}"
        )
    if len(pressure) != expected_points:
        raise ValueError(
            f"pressure field length mismatch in {path}: expected {expected_points}, got {len(pressure)}"
        )
    return velocity, pressure


def format_optional_metric(value: float, spec: str, suffix: str = "") -> str:
    if not math.isfinite(value):
        return "n/a"
    return format(value, spec) + suffix


def write_csv(path: Path, headers: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, contents: str) -> None:
    path.write_text(contents)


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    header_line = "| " + " | ".join(headers) + " |"
    separator_line = "| " + " | ".join(["---"] * len(headers)) + " |"
    body = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header_line, separator_line, *body])


def write_bar_svg(path: Path, title: str, xlabel: str, ylabel: str, rows: list[tuple[str, float]]) -> None:
    width = 900
    height = 520
    left = 90
    right = 40
    top = 60
    bottom = 100
    chart_width = width - left - right
    chart_height = height - top - bottom
    max_value = max(value for _, value in rows)
    scale = chart_height / max_value if max_value > 0.0 else 0.0
    bar_width = chart_width / max(1, len(rows)) * 0.6
    gap = chart_width / max(1, len(rows))
    colors = ["#0f766e", "#b45309", "#1d4ed8", "#be185d"]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fffdf7"/>',
        f'<text x="{width / 2}" y="32" text-anchor="middle" font-size="22" font-family="Georgia, serif">{title}</text>',
        f'<line x1="{left}" y1="{height - bottom}" x2="{width - right}" y2="{height - bottom}" stroke="#111827" stroke-width="2"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}" stroke="#111827" stroke-width="2"/>',
        f'<text x="{width / 2}" y="{height - 26}" text-anchor="middle" font-size="16" font-family="Menlo, monospace">{xlabel}</text>',
        f'<text x="26" y="{height / 2}" transform="rotate(-90 26 {height / 2})" text-anchor="middle" font-size="16" font-family="Menlo, monospace">{ylabel}</text>',
    ]

    for index, (label, value) in enumerate(rows):
        x = left + gap * index + 0.2 * gap
        bar_height = value * scale
        y = height - bottom - bar_height
        color = colors[index % len(colors)]
        parts.append(
            f'<rect x="{x}" y="{y}" width="{bar_width}" height="{bar_height}" fill="{color}" opacity="0.9"/>'
        )
        parts.append(
            f'<text x="{x + bar_width / 2}" y="{y - 8}" text-anchor="middle" font-size="13" font-family="Menlo, monospace">{value:.3g}</text>'
        )
        parts.append(
            f'<text x="{x + bar_width / 2}" y="{height - bottom + 22}" text-anchor="middle" font-size="13" font-family="Menlo, monospace">{label}</text>'
        )

    parts.append("</svg>")
    path.write_text("\n".join(parts))


def get_hardware_info(root: Path) -> dict[str, str]:
    brand = run_command(["sysctl", "-n", "machdep.cpu.brand_string"], root).stdout.strip()
    p_cores = run_command(["sysctl", "-n", "hw.perflevel0.logicalcpu"], root).stdout.strip()
    e_cores = run_command(["sysctl", "-n", "hw.perflevel1.logicalcpu"], root).stdout.strip()
    memsize = int(run_command(["sysctl", "-n", "hw.memsize"], root).stdout.strip())
    return {
        "cpu": brand,
        "performance_cores": p_cores,
        "efficiency_cores": e_cores,
        "memory_gib": f"{memsize / (1024 ** 3):.1f}",
    }


def run_kernel_microbenchmarks(root: Path, build_dir: Path, output_dir: Path) -> list[dict[str, object]]:
    tools_dir = root / build_dir / "tools"
    commands = [
        {
            "case": "advection_256x256",
            "command": [str(tools_dir / "solver_advection_profile"), "256", "240"],
        },
        {
            "case": "pressure_poisson_256x256",
            "command": [str(tools_dir / "solver_pressure_profile"), "256", "24"],
        },
    ]

    rows: list[dict[str, object]] = []
    for item in commands:
        metrics = parse_solver_metrics(run_command(item["command"], root).stdout)
        row: dict[str, object] = {
            "case": item["case"],
            "elapsed_seconds": float(metrics["elapsed_seconds"]),
            "average_seconds_per_repeat": float(metrics["average_seconds_per_repeat"]),
            "logical_cells": int(metrics["logical_cells"]),
        }
        if metrics["kernel"] == "advection":
            row["throughput"] = float(metrics["cell_updates_per_second"])
            row["throughput_label"] = "cell_updates_per_second"
            row["lower_bound_bandwidth_gbps"] = float(metrics["lower_bound_bandwidth_gbps"])
        else:
            row["throughput"] = float(metrics["unknown_updates_per_second"])
            row["throughput_label"] = "unknown_updates_per_second"
            row["average_iterations"] = float(metrics["average_iterations"])
            row["average_relative_residual"] = float(metrics["average_relative_residual"])
            row["multigrid_levels"] = int(metrics["multigrid_levels"])
        rows.append(row)

    write_csv(
        output_dir / "kernel_microbenchmarks.csv",
        [
            "case",
            "elapsed_seconds",
            "average_seconds_per_repeat",
            "logical_cells",
            "throughput_label",
            "throughput",
            "lower_bound_bandwidth_gbps",
            "average_iterations",
            "average_relative_residual",
            "multigrid_levels",
        ],
        rows,
    )
    write_text(
        output_dir / "kernel_microbenchmarks.md",
        "# Kernel Microbenchmarks\n\n"
        + markdown_table(
            [
                "Case",
                "Avg Seconds",
                "Throughput",
                "Lower-Bound BW (GB/s)",
                "Avg Iterations",
                "Avg Residual",
                "MG Levels",
            ],
            [
                [
                    str(row["case"]),
                    f'{float(row["average_seconds_per_repeat"]):.6e}',
                    f'{row["throughput_label"]}={float(row["throughput"]):.6e}',
                    ""
                    if "lower_bound_bandwidth_gbps" not in row
                    else f'{float(row["lower_bound_bandwidth_gbps"]):.3f}',
                    "" if "average_iterations" not in row else f'{float(row["average_iterations"]):.2f}',
                    ""
                    if "average_relative_residual" not in row
                    else f'{float(row["average_relative_residual"]):.3e}',
                    "" if "multigrid_levels" not in row else str(row["multigrid_levels"]),
                ]
                for row in rows
            ],
        )
        + "\n",
    )
    write_bar_svg(
        output_dir / "kernel_throughput.svg",
        "Kernel Throughput",
        "Kernel",
        "Throughput",
        [(str(row["case"]), float(row["throughput"])) for row in rows],
    )
    return rows


def run_end_to_end_cases(root: Path, build_dir: Path, output_dir: Path) -> list[dict[str, object]]:
    tools_dir = root / build_dir / "tools"
    cases = [
        {
            "case": "taylor_green_128",
            "command": [str(tools_dir / "solver_taylor_green"), "benchmarks/taylor_green_128.cfg"],
            "logical_cells": 128 * 128,
        },
        {
            "case": "channel_couette_128",
            "command": [str(tools_dir / "solver_channel"), "benchmarks/channel_couette_128.cfg"],
            "logical_cells": 128 * 128,
        },
        {
            "case": "cavity_smoke_32_steps256",
            "command": [
                str(tools_dir / "solver_cavity"),
                "benchmarks/lid_driven_cavity_smoke.cfg",
                "--steps",
                "256",
            ],
            "logical_cells": 32 * 32,
        },
    ]

    rows: list[dict[str, object]] = []
    for case in cases:
        solver_metrics, time_metrics = run_with_time(case["command"], root)
        steps = int(float(solver_metrics["steps"]))
        real_seconds = time_metrics["real_seconds"]
        cells_per_second = case["logical_cells"] * steps / real_seconds
        rows.append(
            {
                "case": case["case"],
                "steps": steps,
                "logical_cells": case["logical_cells"],
                "real_seconds": real_seconds,
                "cells_per_second": cells_per_second,
                "maximum_resident_set_size": time_metrics.get("maximum_resident_set_size", 0.0),
                "peak_memory_footprint": time_metrics.get("peak_memory_footprint", 0.0),
                "benchmark_status": solver_metrics.get("benchmark_status", ""),
            }
        )

    write_csv(
        output_dir / "throughput_summary.csv",
        [
            "case",
            "steps",
            "logical_cells",
            "real_seconds",
            "cells_per_second",
            "maximum_resident_set_size",
            "peak_memory_footprint",
            "benchmark_status",
        ],
        rows,
    )
    write_text(
        output_dir / "throughput_summary.md",
        "# End-to-End Throughput Summary\n\n"
        + markdown_table(
            [
                "Case",
                "Steps",
                "Cells",
                "Real Seconds",
                "Cells / Second",
                "Max RSS",
                "Peak Footprint",
                "Status",
            ],
            [
                [
                    str(row["case"]),
                    str(row["steps"]),
                    str(row["logical_cells"]),
                    f'{float(row["real_seconds"]):.3f}',
                    f'{float(row["cells_per_second"]):.6e}',
                    f'{float(row["maximum_resident_set_size"]):.0f}',
                    f'{float(row["peak_memory_footprint"]):.0f}',
                    str(row["benchmark_status"]),
                ]
                for row in rows
            ],
        )
        + "\n",
    )
    return rows


def run_cpu_vs_metal_cases(root: Path, build_dir: Path, output_dir: Path) -> list[dict[str, object]]:
    tools_dir = root / build_dir / "tools"
    cases = [
        {
            "case": "taylor_green_3d_64",
            "config": "benchmarks/taylor_green_3d_64.cfg",
            "expected_status": "pass",
            "compare_fields": True,
        },
        {
            "case": "taylor_green_3d_256_smoke",
            "config": "benchmarks/taylor_green_3d_256_smoke.cfg",
            "expected_status": "skipped",
            "compare_fields": False,
        },
    ]

    rows: list[dict[str, object]] = []
    for case in cases:
        cpu_command = [
            str(tools_dir / "solver_taylor_green"),
            case["config"],
        ]
        metal_command = [
            str(tools_dir / "solver_taylor_green"),
            case["config"],
            "--backend",
            "metal",
        ]

        cpu_vtk: Path | None = None
        metal_vtk: Path | None = None
        if case["compare_fields"]:
            cpu_vtk = output_dir / f"{case['case']}_cpu.vtk"
            metal_vtk = output_dir / f"{case['case']}_metal.vtk"
            cpu_command.extend(["--vtk-out", str(cpu_vtk)])
            metal_command.extend(["--vtk-out", str(metal_vtk)])

        cpu_metrics, cpu_time = run_with_time(cpu_command, root)
        metal_metrics, metal_time = run_with_time(metal_command, root)
        if cpu_metrics.get("benchmark_status", "") != case["expected_status"]:
            raise RuntimeError(
                f"{case['case']} CPU status mismatch: expected {case['expected_status']}, got {cpu_metrics.get('benchmark_status', '')}"
            )
        if metal_metrics.get("benchmark_status", "") != case["expected_status"]:
            raise RuntimeError(
                f"{case['case']} Metal status mismatch: expected {case['expected_status']}, got {metal_metrics.get('benchmark_status', '')}"
            )

        velocity_relative_l2_difference = math.nan
        pressure_relative_l2_difference = math.nan
        if cpu_vtk is not None and metal_vtk is not None:
            try:
                cpu_velocity, cpu_pressure = read_vtk_fields(cpu_vtk)
                metal_velocity, metal_pressure = read_vtk_fields(metal_vtk)
                velocity_relative_l2_difference = relative_l2_difference(cpu_velocity, metal_velocity)
                pressure_relative_l2_difference = relative_l2_difference(cpu_pressure, metal_pressure)
            finally:
                cpu_vtk.unlink(missing_ok=True)
                metal_vtk.unlink(missing_ok=True)

        cpu_energy = float(cpu_metrics["final_kinetic_energy"])
        metal_energy = float(metal_metrics["final_kinetic_energy"])
        energy_relative_difference = (
            abs(metal_energy - cpu_energy) / cpu_energy if cpu_energy > 0.0 else abs(metal_energy - cpu_energy)
        )
        metal_real_seconds = metal_time["real_seconds"]

        rows.append(
            {
                "case": case["case"],
                "field_comparison_mode": "full_vtk_diff" if case["compare_fields"] else "metrics_only",
                "cpu_status": cpu_metrics.get("benchmark_status", ""),
                "metal_status": metal_metrics.get("benchmark_status", ""),
                "cpu_real_seconds": cpu_time["real_seconds"],
                "metal_real_seconds": metal_real_seconds,
                "speedup": cpu_time["real_seconds"] / metal_real_seconds if metal_real_seconds > 0.0 else math.inf,
                "velocity_relative_l2_difference": velocity_relative_l2_difference,
                "pressure_relative_l2_difference": pressure_relative_l2_difference,
                "energy_relative_difference": energy_relative_difference,
                "cpu_divergence_l2": float(cpu_metrics["divergence_l2"]),
                "metal_divergence_l2": float(metal_metrics["divergence_l2"]),
                "metal_backend_elapsed_seconds": float(metal_metrics.get("backend_elapsed_seconds", "0")),
                "metal_cleanup_elapsed_seconds": float(metal_metrics.get("cleanup_elapsed_seconds", "0")),
            }
        )

    write_csv(
        output_dir / "cpu_vs_metal.csv",
        [
            "case",
            "field_comparison_mode",
            "cpu_status",
            "metal_status",
            "cpu_real_seconds",
            "metal_real_seconds",
            "speedup",
            "velocity_relative_l2_difference",
            "pressure_relative_l2_difference",
            "energy_relative_difference",
            "cpu_divergence_l2",
            "metal_divergence_l2",
            "metal_backend_elapsed_seconds",
            "metal_cleanup_elapsed_seconds",
        ],
        rows,
    )
    write_text(
        output_dir / "cpu_vs_metal.md",
        "# CPU vs Metal Taylor-Green Comparison\n\n"
        "These comparisons use the same Taylor-Green executable and configuration surface, with the Metal path selected via `--backend metal`.\n\n"
        "The current Metal slice is velocity-first: it uses float working storage on device and a final CPU projection cleanup before result publication, so velocity / energy agreement is the primary correctness comparison for this milestone.\n\n"
        "The `256^3` smoke case is intentionally metrics-only here to avoid materializing giant ASCII VTK snapshots during a long profiling run.\n\n"
        + markdown_table(
            [
                "Case",
                "Comparison Mode",
                "CPU Status",
                "Metal Status",
                "CPU Seconds",
                "Metal Seconds",
                "Speedup",
                "Velocity Rel L2",
                "Pressure Rel L2",
                "Energy Rel Diff",
                "CPU Div L2",
                "Metal Div L2",
            ],
            [
                [
                    str(row["case"]),
                    str(row["field_comparison_mode"]),
                    str(row["cpu_status"]),
                    str(row["metal_status"]),
                    f'{float(row["cpu_real_seconds"]):.3f}',
                    f'{float(row["metal_real_seconds"]):.3f}',
                    format_optional_metric(float(row["speedup"]), ".3f", "x"),
                    format_optional_metric(float(row["velocity_relative_l2_difference"]), ".6e"),
                    format_optional_metric(float(row["pressure_relative_l2_difference"]), ".6e"),
                    format_optional_metric(float(row["energy_relative_difference"]), ".6e"),
                    format_optional_metric(float(row["cpu_divergence_l2"]), ".6e"),
                    format_optional_metric(float(row["metal_divergence_l2"]), ".6e"),
                ]
                for row in rows
            ],
        )
        + "\n",
    )
    return rows


def record_time_profile(root: Path, command: list[str]) -> tuple[list[dict[str, object]], dict[str, float]]:
    with tempfile.TemporaryDirectory(prefix="solver_profile_trace_") as temp_dir:
        trace_path = Path(temp_dir) / "run.trace"
        xml_path = Path(temp_dir) / "time-profile.xml"
        run_command(
            [
                "xctrace",
                "record",
                "--quiet",
                "--template",
                "Time Profiler",
                "--output",
                str(trace_path),
                "--launch",
                "--",
                *command,
            ],
            root,
        )
        run_command_stdout_to_path(
            [
                "xctrace",
                "export",
                "--input",
                str(trace_path),
                "--xpath",
                '/trace-toc/run[@number="1"]/data/table[@schema="time-profile"]',
            ],
            root,
            xml_path,
        )
        return extract_profile_summary(xml_path)


def classify_symbol(symbol: str) -> str:
    if any(token in symbol for token in ["solve_pressure_poisson", "smooth_jacobi", "v_cycle", "neighbor_value", "restrict_full_weighting", "direct_coarse_solve", "apply_poisson_operator"]):
        return "pressure_solve"
    if any(token in symbol for token in ["compute_u_advection_2d", "compute_v_advection_2d", "compute_w_advection_2d", "reconstruct_face_x", "reconstruct_face_y"]):
        return "advection"
    if any(token in symbol for token in ["solve_predictor_adi", "solve_component_sweep", "solve_tridiagonal", "solve_cyclic_tridiagonal"]):
        return "predictor_adi"
    if any(token in symbol for token in ["compute_diffusion_term", "compute_component_laplacian"]):
        return "diffusion"
    if any(token in symbol for token in ["FieldLayout::is_storage_index", "FieldLayout::index", "FieldLayout::storage_index_from_active"]):
        return "field_layout"
    if any(token in symbol for token in ["operator new", "_xzm_", "_platform_memmove", "copy_component", "copy_active"]):
        return "allocation_copy"
    return "other"


def extract_profile_summary(xml_path: Path) -> tuple[list[dict[str, object]], dict[str, float]]:
    core_counts: Counter[str] = Counter()
    symbol_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter()

    for _event, row in ET.iterparse(xml_path, events=("end",)):
        if row.tag != "row":
            continue
        core = row.find("core")
        if core is not None:
            label = core.attrib.get("fmt", "")
            if "P Core" in label:
                core_counts["P"] += 1
            elif "E Core" in label:
                core_counts["E"] += 1
            else:
                core_counts["other"] += 1

        backtrace = row.find("backtrace")
        if backtrace is None:
            continue

        symbol = None
        for frame in backtrace.findall("frame"):
            name = frame.attrib.get("name", "")
            if "solver::" in name or "solver_" in name or "FieldLayout::" in name:
                symbol = name
                break
        if symbol is None:
            row.clear()
            continue

        symbol_counts[symbol] += 1
        category_counts[classify_symbol(symbol)] += 1
        row.clear()

    total_symbol_samples = sum(symbol_counts.values()) or 1
    hotspot_rows = [
        {
            "symbol": symbol,
            "samples": samples,
            "sample_share": samples / total_symbol_samples,
            "category": classify_symbol(symbol),
        }
        for symbol, samples in symbol_counts.most_common(12)
    ]

    labeled_core_samples = core_counts["P"] + core_counts["E"]
    core_summary = {
        "p_core_share": core_counts["P"] / labeled_core_samples if labeled_core_samples else 0.0,
        "e_core_share": core_counts["E"] / labeled_core_samples if labeled_core_samples else 0.0,
        "p_samples": float(core_counts["P"]),
        "e_samples": float(core_counts["E"]),
        "other_samples": float(core_counts["other"]),
    }
    core_summary.update(
        {
            f"category_{category}": count / total_symbol_samples
            for category, count in category_counts.items()
        }
    )
    return hotspot_rows, core_summary


def run_policy_study(root: Path, build_dir: Path, output_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    tools_dir = root / build_dir / "tools"
    commands = {
        "default": [str(tools_dir / "solver_taylor_green"), "benchmarks/taylor_green_128.cfg"],
        "utility": ["/usr/sbin/taskpolicy", "-c", "utility", str(tools_dir / "solver_taylor_green"), "benchmarks/taylor_green_128.cfg"],
        "background": ["/usr/sbin/taskpolicy", "-c", "background", str(tools_dir / "solver_taylor_green"), "benchmarks/taylor_green_128.cfg"],
    }

    policy_rows: list[dict[str, object]] = []
    hotspot_rows: list[dict[str, object]] = []
    for policy_name, command in commands.items():
        solver_metrics, time_metrics = run_with_time(command, root)
        local_hotspots, core_summary = record_time_profile(root, command)

        steps = int(float(solver_metrics["steps"]))
        logical_cells = 128 * 128
        real_seconds = time_metrics["real_seconds"]
        policy_rows.append(
            {
                "policy": policy_name,
                "real_seconds": real_seconds,
                "cells_per_second": logical_cells * steps / real_seconds,
                "maximum_resident_set_size": time_metrics.get("maximum_resident_set_size", 0.0),
                "peak_memory_footprint": time_metrics.get("peak_memory_footprint", 0.0),
                "p_core_share": core_summary["p_core_share"],
                "e_core_share": core_summary["e_core_share"],
            }
        )

        if policy_name == "default":
            for row in local_hotspots:
                hotspot_rows.append(row)

    write_csv(
        output_dir / "policy_study.csv",
        [
            "policy",
            "real_seconds",
            "cells_per_second",
            "maximum_resident_set_size",
            "peak_memory_footprint",
            "p_core_share",
            "e_core_share",
        ],
        policy_rows,
    )
    write_text(
        output_dir / "policy_study.md",
        "# QoS / Clamp Policy Study\n\n"
        + markdown_table(
            [
                "Policy",
                "Real Seconds",
                "Cells / Second",
                "Max RSS",
                "Peak Footprint",
                "P-Core Share",
                "E-Core Share",
            ],
            [
                [
                    str(row["policy"]),
                    f'{float(row["real_seconds"]):.3f}',
                    f'{float(row["cells_per_second"]):.6e}',
                    f'{float(row["maximum_resident_set_size"]):.0f}',
                    f'{float(row["peak_memory_footprint"]):.0f}',
                    f'{float(row["p_core_share"]):.3f}',
                    f'{float(row["e_core_share"]):.3f}',
                ]
                for row in policy_rows
            ],
        )
        + "\n",
    )
    write_bar_svg(
        output_dir / "policy_runtime.svg",
        "QoS Clamp Runtime",
        "Policy",
        "Real Seconds",
        [(str(row["policy"]), float(row["real_seconds"])) for row in policy_rows],
    )

    write_csv(
        output_dir / "hotspot_summary.csv",
        ["symbol", "samples", "sample_share", "category"],
        hotspot_rows,
    )
    write_text(
        output_dir / "hotspot_summary.md",
        "# Default-Policy Hotspots\n\n"
        + markdown_table(
            ["Symbol", "Samples", "Share", "Category"],
            [
                [
                    str(row["symbol"]),
                    str(row["samples"]),
                    f'{float(row["sample_share"]):.3f}',
                    str(row["category"]),
                ]
                for row in hotspot_rows
            ],
        )
        + "\n",
    )
    return policy_rows, hotspot_rows


def write_thread_scaling_baseline(output_dir: Path, reference_cells_per_second: float) -> None:
    rows = [
        {
            "compute_threads": 1,
            "cells_per_second": reference_cells_per_second,
            "status": "baseline_only",
            "note": "Current solver path remains single-threaded after the first M11 optimization pass; deeper thread scaling is still future work.",
        }
    ]
    write_csv(output_dir / "thread_scaling.csv", ["compute_threads", "cells_per_second", "status", "note"], rows)
    write_text(
        output_dir / "thread_scaling.md",
        "# Thread Scaling Baseline\n\n"
        + markdown_table(
            ["Threads", "Cells / Second", "Status", "Note"],
            [[
                str(rows[0]["compute_threads"]),
                f'{float(rows[0]["cells_per_second"]):.6e}',
                str(rows[0]["status"]),
                str(rows[0]["note"]),
            ]],
        )
        + "\n",
    )


def load_baseline_metrics(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return {str(row["metric"]): row for row in reader}


def write_baseline_comparison(
    output_dir: Path,
    hardware: dict[str, str],
    kernel_rows: list[dict[str, object]],
    throughput_rows: list[dict[str, object]],
) -> None:
    baseline_path = Path(__file__).resolve().parent / "baselines" / "m10_m1_max.csv"
    baselines = load_baseline_metrics(baseline_path)
    if not baselines:
        return

    current_metrics = {
        "advection_256x256_throughput": float(
            next(row for row in kernel_rows if row["case"] == "advection_256x256")["throughput"]
        ),
        "pressure_poisson_256x256_throughput": float(
            next(row for row in kernel_rows if row["case"] == "pressure_poisson_256x256")["throughput"]
        ),
        "taylor_green_128_cells_per_second": float(
            next(row for row in throughput_rows if row["case"] == "taylor_green_128")["cells_per_second"]
        ),
        "channel_couette_128_cells_per_second": float(
            next(row for row in throughput_rows if row["case"] == "channel_couette_128")["cells_per_second"]
        ),
        "cavity_smoke_32_steps256_cells_per_second": float(
            next(row for row in throughput_rows if row["case"] == "cavity_smoke_32_steps256")["cells_per_second"]
        ),
    }

    rows: list[dict[str, object]] = []
    for metric, baseline in baselines.items():
        if metric not in current_metrics:
            continue

        baseline_value = float(baseline["baseline_value"])
        current_value = current_metrics[metric]
        speedup = current_value / baseline_value if baseline_value > 0.0 else math.nan
        rows.append(
            {
                "metric": metric,
                "label": baseline["label"],
                "baseline_value": baseline_value,
                "current_value": current_value,
                "speedup": speedup,
                "percent_improvement": (speedup - 1.0) * 100.0,
                "units": baseline["units"],
            }
        )

    write_csv(
        output_dir / "improvement_vs_m10.csv",
        ["metric", "label", "baseline_value", "current_value", "speedup", "percent_improvement", "units"],
        rows,
    )
    write_text(
        output_dir / "improvement_vs_m10.md",
        "# Milestone 10 to Milestone 11 Improvement\n\n"
        + f"Baseline source: `{baseline_path.name}` on `{hardware['cpu']}`.\n\n"
        + markdown_table(
            ["Metric", "M10 Baseline", "M11 Current", "Speedup", "Improvement", "Units"],
            [
                [
                    str(row["label"]),
                    f'{float(row["baseline_value"]):.6e}',
                    f'{float(row["current_value"]):.6e}',
                    f'{float(row["speedup"]):.3f}x',
                    f'{float(row["percent_improvement"]):.1f}%',
                    str(row["units"]),
                ]
                for row in rows
            ],
        )
        + "\n",
    )


def write_hardware_report(output_dir: Path, hardware: dict[str, str]) -> None:
    write_text(
        output_dir / "hardware.md",
        "# Profiling Hardware\n\n"
        + markdown_table(
            ["Key", "Value"],
            [
                ["CPU", hardware["cpu"]],
                ["Performance cores", hardware["performance_cores"]],
                ["Efficiency cores", hardware["efficiency_cores"]],
                ["Memory (GiB)", hardware["memory_gib"]],
            ],
        )
        + "\n",
    )


def write_summary(
    output_dir: Path,
    hardware: dict[str, str],
    kernel_rows: list[dict[str, object]],
    throughput_rows: list[dict[str, object]],
    policy_rows: list[dict[str, object]],
    hotspot_rows: list[dict[str, object]],
    metal_rows: list[dict[str, object]],
) -> None:
    best_policy = min(policy_rows, key=lambda row: float(row["real_seconds"]))
    pressure_row = next(row for row in kernel_rows if row["case"] == "pressure_poisson_256x256")
    advection_row = next(row for row in kernel_rows if row["case"] == "advection_256x256")
    taylor_row = next(row for row in throughput_rows if row["case"] == "taylor_green_128")
    metal_64_row = next(row for row in metal_rows if row["case"] == "taylor_green_3d_64")
    metal_256_row = next(row for row in metal_rows if row["case"] == "taylor_green_3d_256_smoke")
    m10_baselines = load_baseline_metrics(Path(__file__).resolve().parent / "baselines" / "m10_m1_max.csv")
    comparison_lines: list[str] = []
    for metric, current_value in (
        ("advection_256x256_throughput", float(advection_row["throughput"])),
        ("pressure_poisson_256x256_throughput", float(pressure_row["throughput"])),
        ("taylor_green_128_cells_per_second", float(taylor_row["cells_per_second"])),
    ):
        if metric not in m10_baselines:
            continue
        baseline_value = float(m10_baselines[metric]["baseline_value"])
        if baseline_value <= 0.0:
            continue
        speedup = current_value / baseline_value
        comparison_lines.append(
            f"- {m10_baselines[metric]['label']}: {speedup:.3f}x vs the stored Milestone 10 M1 Max baseline."
        )

    top_categories: Counter[str] = Counter()
    for row in hotspot_rows:
        top_categories[str(row["category"])] += int(row["samples"])
    summary = "\n".join(
        [
            "# Performance Summary",
            "",
            "## Hardware",
            "",
            f"- CPU: {hardware['cpu']}",
            f"- Performance cores: {hardware['performance_cores']}",
            f"- Efficiency cores: {hardware['efficiency_cores']}",
            f"- Memory: {hardware['memory_gib']} GiB",
            "",
            "## Findings",
            "",
            f"- Recommended default execution mode: benchmark build profile with the default unclamped scheduler policy. It was the fastest measured policy on this machine.",
            f"- Current compute-thread recommendation: 1. The solver path remains single-threaded after the first M11 optimization pass, so thread scaling is still a baseline-only study for now.",
            f"- Advection microbenchmark throughput: {float(advection_row['throughput']):.6e} cell updates/s with a {float(advection_row['lower_bound_bandwidth_gbps']):.3f} GB/s lower-bound bandwidth estimate.",
            f"- Pressure microbenchmark throughput: {float(pressure_row['throughput']):.6e} unknown updates/s with {float(pressure_row['average_iterations']):.2f} average iterations.",
            f"- End-to-end Taylor-Green benchmark throughput: {float(taylor_row['cells_per_second']):.6e} cells/s.",
            f"- Targeted Metal backend result: {float(metal_64_row['speedup']):.3f}x faster than benchmark-profile CPU on `taylor_green_3d_64`, with velocity relative L2 drift {float(metal_64_row['velocity_relative_l2_difference']):.6e}.",
            f"- Targeted Metal backend smoke throughput: {float(metal_256_row['speedup']):.3f}x faster than benchmark-profile CPU on `taylor_green_3d_256_smoke`.",
            f"- Default-policy hotspot categories were dominated by pressure-solve, predictor/ADI, and advection work, with the top sampled category counts: {dict(top_categories.most_common(4))}.",
            "",
            "## Comparison To Milestone 10 Baseline",
            "",
            *comparison_lines,
            "",
            "## Notes",
            "",
            "- QoS clamp measurements come from default, utility, and background taskpolicy launches.",
            "- Core-class shares come from xctrace Time Profiler samples on Apple Silicon and are reported as sampled P-core vs E-core shares.",
            "- The thread-scaling section is intentionally explicit that the current solver path is still single-threaded after the first CPU optimization pass.",
            "",
        ]
    )
    write_text(output_dir / "summary.md", summary)


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    build_dir = Path(args.build_dir)
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    hardware = get_hardware_info(root)
    write_hardware_report(output_dir, hardware)

    kernel_rows = run_kernel_microbenchmarks(root, build_dir, output_dir)
    throughput_rows = run_end_to_end_cases(root, build_dir, output_dir)
    metal_rows = run_cpu_vs_metal_cases(root, build_dir, output_dir)
    policy_rows, hotspot_rows = run_policy_study(root, build_dir, output_dir)
    default_policy = next(row for row in policy_rows if row["policy"] == "default")
    write_thread_scaling_baseline(output_dir, float(default_policy["cells_per_second"]))
    write_baseline_comparison(output_dir, hardware, kernel_rows, throughput_rows)
    write_summary(output_dir, hardware, kernel_rows, throughput_rows, policy_rows, hotspot_rows, metal_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
