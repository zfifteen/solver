# Deterministic Benchmark Suite

This file documents the deterministic benchmark suite that Milestone 14 treats
as the manual regression gate. This manual suite is the Milestone 14
benchmark-threshold gate; it is intentionally not enforced by GitHub Actions.

The source of truth is
[`benchmark_suite.csv`](./benchmark_suite.csv).
The validation harness in
[`run_validation_suite.py`](./run_validation_suite.py)
reads that manifest directly.

GitHub Actions does not run this suite. The workflow in
[`deterministic-validation.yml`](../.github/workflows/deterministic-validation.yml)
is intentionally limited to fast deterministic build-and-test checks, while the
full validation harness remains a manual developer workflow:

```bash
./validation/run_validation_suite.py --build-dir build/deterministic --output-dir validation/latest
```

The current suite covers:

- Couette validation at `128 x 128`
- Poiseuille validation at `128 x 128`
- Re = 100 lid-driven cavity validation at `128 x 128`
- 2D Taylor-Green decay validation at `128 x 128`
- 3D Taylor-Green decay validation at `64 x 64 x 64`

Each suite entry records:

- the executable to run
- the benchmark configuration file
- the primary reported metric
- the acceptance threshold text used in the generated report bundle
