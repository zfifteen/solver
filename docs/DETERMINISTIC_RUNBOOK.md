# Deterministic Runbook

This runbook defines the raw command sequence for the Milestone 15 release-candidate workflow. It is the operator-facing companion to the automated wrapper in [`tools/run_release_candidate_suite.sh`](../tools/run_release_candidate_suite.sh).

## Supported Scope

- Platform: macOS on Apple Silicon
- Deterministic reference build: `deterministic`
- Performance/profiling build: `benchmark`
- Restartable full-simulation path: lid-driven cavity
- Supported Metal scope: 3D periodic Taylor-Green only

## Deterministic Configure / Build / Test

Configure the deterministic build:

```bash
cmake --preset deterministic
```

Build the deterministic targets:

```bash
cmake --build --preset deterministic -j 3
```

Run the deterministic test suite:

```bash
ctest --preset deterministic
```

## Full Deterministic Validation

Refresh the canonical validation artifacts:

```bash
./validation/run_validation_suite.py --build-dir build/deterministic --output-dir validation/latest
```

The primary acceptance summaries are:

- [`validation/latest/summary.md`](../validation/latest/summary.md)
- [`validation/latest/benchmark_summary.md`](../validation/latest/benchmark_summary.md)

## Checkpoint / Restart Verification

Write a short cavity checkpoint:

```bash
build/deterministic/tools/solver_cavity benchmarks/lid_driven_cavity_smoke.cfg --steps 6 --checkpoint-out /tmp/solver.chk
```

Restart from that checkpoint:

```bash
build/deterministic/tools/solver_cavity benchmarks/lid_driven_cavity_smoke.cfg --restart /tmp/solver.chk --steps 6
```

This is the supported manual restart smoke for the current release candidate. The automated runtime tests remain the main compatibility proof for the checkpoint format.

## Benchmark Configure / Build

Configure the benchmark build:

```bash
cmake --preset benchmark
```

Build the benchmark targets:

```bash
cmake --build --preset benchmark -j 3
```

## Profiling Report Generation

Refresh the canonical profiling artifacts:

```bash
./profiling/run_profile_suite.py --build-dir build/benchmark --output-dir profiling/latest
```

The primary performance-report summaries are:

- [`profiling/latest/summary.md`](../profiling/latest/summary.md)
- [`profiling/latest/throughput_summary.md`](../profiling/latest/throughput_summary.md)
- [`profiling/latest/hotspot_summary.md`](../profiling/latest/hotspot_summary.md)

If you need to detach the long `xctrace` pass from an interactive terminal, the convenience helper remains available:

```bash
tools/profile_suite_background.sh build/benchmark profiling/latest
```

That helper is optional. The canonical release-candidate evidence path is the foreground profiling harness command above.
