#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the deterministic validation suite.")
    parser.add_argument(
        "--build-dir",
        default="build/deterministic",
        help="Build directory that contains the compiled solver tools.",
    )
    parser.add_argument(
        "--output-dir",
        default="validation/latest",
        help="Directory where reports, tables, and plots should be written.",
    )
    return parser.parse_args()


def run_command(args: list[str], cwd: Path) -> str:
    completed = subprocess.run(
        args,
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout


def parse_metrics(stdout: str) -> dict[str, str]:
    metrics: dict[str, str] = {}
    for line in stdout.splitlines():
        if ": " not in line:
            continue
        key, value = line.split(": ", 1)
        metrics[key.strip()] = value.strip()
    return metrics


def load_benchmark_suite(root: Path) -> list[dict[str, str]]:
    manifest_path = root / "validation/benchmark_suite.csv"
    with manifest_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(metrics: dict[str, str], key: str) -> float:
    return float(metrics[key])


def observed_order(coarse_error: float, fine_error: float) -> float:
    return math.log(coarse_error / fine_error, 2.0)


def write_csv(path: Path, headers: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    header_line = "| " + " | ".join(headers) + " |"
    separator_line = "| " + " | ".join(["---"] * len(headers)) + " |"
    body = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header_line, separator_line, *body])


def write_text(path: Path, contents: str) -> None:
    path.write_text(contents)


def write_loglog_svg(
    path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    series: list[tuple[str, list[tuple[float, float]]]],
) -> None:
    width = 900
    height = 560
    left = 100
    right = 40
    top = 60
    bottom = 80
    colors = ["#0f766e", "#b45309", "#1d4ed8", "#be185d"]

    xs = [math.log10(x) for _, points in series for x, _ in points if x > 0.0]
    ys = [math.log10(y) for _, points in series for _, y in points if y > 0.0]
    if not xs or not ys:
        raise ValueError("write_loglog_svg needs strictly positive data")

    x_min = min(xs)
    x_max = max(xs)
    y_min = min(ys)
    y_max = max(ys)
    if x_min == x_max:
        x_min -= 0.5
        x_max += 0.5
    if y_min == y_max:
        y_min -= 0.5
        y_max += 0.5

    def map_x(value: float) -> float:
        return left + (math.log10(value) - x_min) * (width - left - right) / (x_max - x_min)

    def map_y(value: float) -> float:
        return height - bottom - (math.log10(value) - y_min) * (height - top - bottom) / (y_max - y_min)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fffdf7"/>',
        f'<text x="{width / 2}" y="32" text-anchor="middle" font-size="22" font-family="Georgia, serif">{title}</text>',
        f'<line x1="{left}" y1="{height - bottom}" x2="{width - right}" y2="{height - bottom}" stroke="#111827" stroke-width="2"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}" stroke="#111827" stroke-width="2"/>',
        f'<text x="{width / 2}" y="{height - 22}" text-anchor="middle" font-size="16" font-family="Menlo, monospace">{xlabel}</text>',
        f'<text x="26" y="{height / 2}" transform="rotate(-90 26 {height / 2})" text-anchor="middle" font-size="16" font-family="Menlo, monospace">{ylabel}</text>',
    ]

    x_ticks = sorted({x for _, points in series for x, _ in points if x > 0.0})
    y_ticks = sorted({y for _, points in series for _, y in points if y > 0.0})
    for value in x_ticks:
        x = map_x(value)
        parts.append(f'<line x1="{x}" y1="{height - bottom}" x2="{x}" y2="{height - bottom + 6}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{x}" y="{height - bottom + 24}" text-anchor="middle" font-size="13" font-family="Menlo, monospace">{value:.3g}</text>')
    for value in y_ticks:
        y = map_y(value)
        parts.append(f'<line x1="{left - 6}" y1="{y}" x2="{left}" y2="{y}" stroke="#111827" stroke-width="1"/>')
        parts.append(f'<text x="{left - 10}" y="{y + 4}" text-anchor="end" font-size="13" font-family="Menlo, monospace">{value:.3g}</text>')

    legend_y = top + 12
    for index, (label, points) in enumerate(series):
        color = colors[index % len(colors)]
        coords = " ".join(f"{map_x(x)},{map_y(y)}" for x, y in points)
        parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="3" points="{coords}"/>')
        for x, y in points:
            parts.append(f'<circle cx="{map_x(x)}" cy="{map_y(y)}" r="4" fill="{color}"/>')
        legend_x = width - right - 220
        parts.append(f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 24}" y2="{legend_y}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="{legend_x + 32}" y="{legend_y + 4}" font-size="14" font-family="Menlo, monospace">{label}</text>')
        legend_y += 24

    parts.append("</svg>")
    path.write_text("\n".join(parts))


def read_vtk_velocity(path: Path, expected_value_count: int) -> list[float]:
    values: list[float] = []
    reading = False
    for line in path.read_text().splitlines():
        if line.startswith("VECTORS velocity double"):
            reading = True
            continue
        if reading and line.startswith("SCALARS pressure double 1"):
            break
        if reading:
            parts = line.split()
            if len(parts) == 3:
                values.extend(float(part) for part in parts)

    if expected_value_count <= 0:
        raise ValueError("read_vtk_velocity expects a positive expected_value_count")
    if len(values) != expected_value_count:
        raise ValueError(
            f"read_vtk_velocity expected {expected_value_count} values in {path}, "
            f"but parsed {len(values)}"
        )
    return values


def relative_l2_difference(left: list[float], right: list[float]) -> float:
    if not left or not right:
        raise ValueError("relative_l2_difference expects non-empty vectors")
    if len(left) != len(right):
        raise ValueError(
            f"relative_l2_difference expects equal-length vectors, got {len(left)} and {len(right)}"
        )
    numerator = math.sqrt(sum((a - b) * (a - b) for a, b in zip(left, right)) / len(left))
    denominator = math.sqrt(sum(b * b for b in right) / len(right))
    return numerator / denominator


def write_operator_reports(
    root: Path, build_dir: Path, output_dir: Path
) -> tuple[dict[str, float], list[dict[str, object]]]:
    resolutions = [16, 32, 64]
    rows: list[dict[str, object]] = []
    results: list[dict[str, float]] = []
    tool = root / build_dir / "tools/solver_operator_mms"

    for resolution in resolutions:
        metrics = parse_metrics(run_command([str(tool), str(resolution)], root))
        result = {
            "resolution": float(metrics["resolution"]),
            "gradient_l2_error": as_float(metrics, "gradient_l2_error"),
            "divergence_l2_error": as_float(metrics, "divergence_l2_error"),
            "laplacian_l2_error": as_float(metrics, "laplacian_l2_error"),
        }
        results.append(result)

    for index, result in enumerate(results):
        row = {
            "resolution": int(result["resolution"]),
            "gradient_l2_error": result["gradient_l2_error"],
            "divergence_l2_error": result["divergence_l2_error"],
            "laplacian_l2_error": result["laplacian_l2_error"],
            "gradient_order": "",
            "divergence_order": "",
            "laplacian_order": "",
        }
        if index > 0:
            coarse = results[index - 1]
            row["gradient_order"] = observed_order(
                coarse["gradient_l2_error"], result["gradient_l2_error"]
            )
            row["divergence_order"] = observed_order(
                coarse["divergence_l2_error"], result["divergence_l2_error"]
            )
            row["laplacian_order"] = observed_order(
                coarse["laplacian_l2_error"], result["laplacian_l2_error"]
            )
        rows.append(row)

    write_csv(
        output_dir / "operator_spatial.csv",
        [
            "resolution",
            "gradient_l2_error",
            "divergence_l2_error",
            "laplacian_l2_error",
            "gradient_order",
            "divergence_order",
            "laplacian_order",
        ],
        rows,
    )

    write_loglog_svg(
        output_dir / "operator_spatial.svg",
        "Operator MMS Spatial Convergence",
        "Resolution",
        "L2 Error",
        [
            ("gradient", [(result["resolution"], result["gradient_l2_error"]) for result in results]),
            ("divergence", [(result["resolution"], result["divergence_l2_error"]) for result in results]),
            ("laplacian", [(result["resolution"], result["laplacian_l2_error"]) for result in results]),
        ],
    )

    md_rows = []
    for row in rows:
        md_rows.append(
            [
                str(row["resolution"]),
                f'{row["gradient_l2_error"]:.6e}',
                f'{row["divergence_l2_error"]:.6e}',
                f'{row["laplacian_l2_error"]:.6e}',
                "" if row["gradient_order"] == "" else f'{float(row["gradient_order"]):.3f}',
                "" if row["divergence_order"] == "" else f'{float(row["divergence_order"]):.3f}',
                "" if row["laplacian_order"] == "" else f'{float(row["laplacian_order"]):.3f}',
            ]
        )
    write_text(
        output_dir / "operator_spatial.md",
        "# Operator MMS Spatial Convergence\n\n"
        + markdown_table(
            [
                "Resolution",
                "Gradient L2",
                "Divergence L2",
                "Laplacian L2",
                "Gradient Order",
                "Divergence Order",
                "Laplacian Order",
            ],
            md_rows,
        )
        + "\n",
    )

    summary = {
        "gradient_order": float(rows[-1]["gradient_order"]),
        "divergence_order": float(rows[-1]["divergence_order"]),
        "laplacian_order": float(rows[-1]["laplacian_order"]),
    }
    return summary, rows


def write_taylor_temporal_reports(
    root: Path, build_dir: Path, output_dir: Path
) -> tuple[dict[str, float], list[dict[str, object]]]:
    tool = root / build_dir / "tools/solver_taylor_green"
    cfl_limits = [1.0, 0.5, 0.25]
    nx = 128
    ny = 128
    nz = 1
    expected_velocity_value_count = 3 * nx * ny * nz
    rows: list[dict[str, object]] = []
    vectors: list[list[float]] = []

    for cfl_limit in cfl_limits:
        config_path = output_dir / f"taylor_green_temporal_cfl_{str(cfl_limit).replace('.', 'p')}.cfg"
        vtk_path = output_dir / f"taylor_green_temporal_cfl_{str(cfl_limit).replace('.', 'p')}.vtk"
        config_path.write_text(
            "\n".join(
                [
                    "# Generated by validation/run_validation_suite.py",
                    f"nx = {nx}",
                    f"ny = {ny}",
                    f"nz = {nz}",
                    "viscosity = 0.01",
                    f"cfl_limit = {cfl_limit}",
                    "final_time = 0.25",
                    "poisson_max_iterations = 200",
                    "poisson_tolerance = 1e-10",
                    "validate_energy = true",
                ]
            )
            + "\n"
        )
        metrics = parse_metrics(run_command([str(tool), str(config_path), "--vtk-out", str(vtk_path)], root))
        rows.append(
            {
                "cfl_limit": cfl_limit,
                "dt": as_float(metrics, "dt"),
                "normalized_energy_error": as_float(metrics, "normalized_energy_error"),
                "velocity_relative_l2_error": as_float(metrics, "velocity_relative_l2_error"),
                "divergence_l2": as_float(metrics, "divergence_l2"),
            }
        )
        vectors.append(read_vtk_velocity(vtk_path, expected_velocity_value_count))

    rows[0]["self_difference_to_finer"] = relative_l2_difference(vectors[0], vectors[1])
    rows[1]["self_difference_to_finer"] = relative_l2_difference(vectors[1], vectors[2])
    rows[2]["self_difference_to_finer"] = ""
    temporal_order = observed_order(rows[0]["self_difference_to_finer"], rows[1]["self_difference_to_finer"])
    rows[1]["self_convergence_order"] = temporal_order

    write_csv(
        output_dir / "taylor_green_temporal.csv",
        [
            "cfl_limit",
            "dt",
            "normalized_energy_error",
            "velocity_relative_l2_error",
            "divergence_l2",
            "self_difference_to_finer",
            "self_convergence_order",
        ],
        rows,
    )

    write_loglog_svg(
        output_dir / "taylor_green_temporal.svg",
        "Taylor-Green Temporal Refinement",
        "dt",
        "Error",
        [
            ("energy", [(row["dt"], row["normalized_energy_error"]) for row in rows]),
            ("velocity", [(row["dt"], row["velocity_relative_l2_error"]) for row in rows]),
            (
                "self-convergence",
                [
                    (rows[0]["dt"], rows[0]["self_difference_to_finer"]),
                    (rows[1]["dt"], rows[1]["self_difference_to_finer"]),
                ],
            ),
        ],
    )

    md_rows = []
    for row in rows:
        md_rows.append(
            [
                f'{row["cfl_limit"]:.2f}',
                f'{row["dt"]:.6e}',
                f'{row["normalized_energy_error"]:.6e}',
                f'{row["velocity_relative_l2_error"]:.6e}',
                f'{row["divergence_l2"]:.6e}',
                "" if row.get("self_difference_to_finer", "") == "" else f'{float(row["self_difference_to_finer"]):.6e}',
                "" if row.get("self_convergence_order", "") == "" else f'{float(row["self_convergence_order"]):.3f}',
            ]
        )
    write_text(
        output_dir / "taylor_green_temporal.md",
        "# Taylor-Green Temporal Refinement\n\n"
        "Temporal order is measured with a self-convergent final-field comparison on the same 128 x 128 grid.\n\n"
        + markdown_table(
            [
                "CFL",
                "dt",
                "Normalized KE Error",
                "Velocity Relative L2",
                "Divergence L2",
                "Self Difference",
                "Order",
            ],
            md_rows,
        )
        + "\n",
    )

    return {"temporal_order": temporal_order}, rows


def write_benchmark_reports(root: Path, build_dir: Path, output_dir: Path) -> list[dict[str, object]]:
    cases = load_benchmark_suite(root)

    rows: list[dict[str, object]] = []
    for case in cases:
        command = [str(root / build_dir / f'tools/{case["tool"]}'), case["config"]]
        metrics = parse_metrics(run_command(command, root))
        rows.append(
            {
                "case": case["case"],
                "reference_dataset": metrics.get("reference_dataset", ""),
                "metric_key": case["metric_key"],
                "metric_value": as_float(metrics, case["metric_key"]),
                "divergence_l2": as_float(metrics, "divergence_l2"),
                "threshold": case["threshold"],
                "benchmark_status": metrics["benchmark_status"],
            }
        )

    write_csv(
        output_dir / "benchmark_summary.csv",
        ["case", "reference_dataset", "metric_key", "metric_value", "divergence_l2", "threshold", "benchmark_status"],
        rows,
    )

    md_rows = []
    for row in rows:
        md_rows.append(
            [
                str(row["case"]),
                str(row["reference_dataset"]),
                f'{row["metric_key"]}={float(row["metric_value"]):.6e}',
                f'{float(row["divergence_l2"]):.6e}',
                str(row["threshold"]),
                str(row["benchmark_status"]),
            ]
        )
    write_text(
        output_dir / "benchmark_summary.md",
        "# Benchmark Summary\n\n"
        + markdown_table(
            ["Case", "Reference Dataset", "Primary Metric", "Divergence L2", "Threshold", "Status"],
            md_rows,
        )
        + "\n",
    )

    return rows


def write_mass_conservation_report(output_dir: Path, benchmark_rows: list[dict[str, object]]) -> None:
    rows = [
        {
            "case": row["case"],
            "divergence_l2": row["divergence_l2"],
            "threshold": "divergence_l2 <= 1e-10",
            "pass": float(row["divergence_l2"]) <= 1.0e-10,
        }
        for row in benchmark_rows
    ]

    write_csv(output_dir / "mass_conservation.csv", ["case", "divergence_l2", "threshold", "pass"], rows)
    write_text(
        output_dir / "mass_conservation.md",
        "# Mass Conservation Summary\n\n"
        + markdown_table(
            ["Case", "Divergence L2", "Threshold", "Pass"],
            [
                [
                    str(row["case"]),
                    f'{float(row["divergence_l2"]):.6e}',
                    str(row["threshold"]),
                    "yes" if bool(row["pass"]) else "no",
                ]
                for row in rows
            ],
        )
        + "\n",
    )


def write_summary(
    output_dir: Path,
    operator_summary: dict[str, float],
    temporal_summary: dict[str, float],
    benchmark_rows: list[dict[str, object]],
) -> bool:
    operator_pass = all(value >= 1.8 for value in operator_summary.values())
    temporal_pass = temporal_summary["temporal_order"] >= 1.8
    benchmark_pass = all(row["benchmark_status"] == "pass" for row in benchmark_rows)
    divergence_pass = all(float(row["divergence_l2"]) <= 1.0e-10 for row in benchmark_rows)
    suite_pass = operator_pass and temporal_pass and benchmark_pass and divergence_pass

    summary = "\n".join(
        [
            "# Validation Summary",
            "",
            f"Overall status: {'PASS' if suite_pass else 'FAIL'}",
            "",
            "## Gates",
            "",
            f"- Operator spatial order >= 1.8: {'pass' if operator_pass else 'fail'}",
            f"- Taylor-Green temporal self-convergence order >= 1.8: {'pass' if temporal_pass else 'fail'}",
            f"- Benchmark thresholds: {'pass' if benchmark_pass else 'fail'}",
            f"- Mass conservation thresholds: {'pass' if divergence_pass else 'fail'}",
            "",
            "## Notes",
            "",
            "- Spatial convergence is reported from the manufactured-solution operator study.",
            "- Temporal convergence is measured from the 2D Taylor-Green final-field self-convergence study on a fixed 128 x 128 grid.",
            "- Benchmark thresholds mirror the technical specification and benchmark harness contracts already in the repo.",
            "",
        ]
    )
    write_text(output_dir / "summary.md", summary)
    return suite_pass


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    build_dir = Path(args.build_dir)
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    operator_summary, _ = write_operator_reports(root, build_dir, output_dir)
    temporal_summary, _ = write_taylor_temporal_reports(root, build_dir, output_dir)
    benchmark_rows = write_benchmark_reports(root, build_dir, output_dir)
    write_mass_conservation_report(output_dir, benchmark_rows)
    suite_pass = write_summary(output_dir, operator_summary, temporal_summary, benchmark_rows)
    return 0 if suite_pass else 1


if __name__ == "__main__":
    sys.exit(main())
