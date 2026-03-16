# Performance Summary

## Hardware

- CPU: Apple M1 Max
- Performance cores: 8
- Efficiency cores: 2
- Memory: 32.0 GiB

## Findings

- Recommended default execution mode: benchmark build profile with the default unclamped scheduler policy. It was the fastest measured policy on this machine.
- Current compute-thread recommendation: 1. The solver path remains single-threaded after the first M11 optimization pass, so thread scaling is still a baseline-only study for now.
- Advection microbenchmark throughput: 1.882870e+07 cell updates/s with a 1.207 GB/s lower-bound bandwidth estimate.
- Pressure microbenchmark throughput: 1.950190e+06 unknown updates/s with 16.00 average iterations.
- End-to-end Taylor-Green benchmark throughput: 1.138549e+06 cells/s.
- Default-policy hotspot categories were dominated by pressure-solve, predictor/ADI, and advection work, with the top sampled category counts: {'other': 31, 'pressure_solve': 26, 'advection': 20, 'predictor_adi': 9}.

## Comparison To Milestone 10 Baseline

- Advection microbenchmark throughput: 1.482x vs the stored Milestone 10 M1 Max baseline.
- Pressure microbenchmark throughput: 6.608x vs the stored Milestone 10 M1 Max baseline.
- Taylor-Green end-to-end throughput: 3.915x vs the stored Milestone 10 M1 Max baseline.

## Notes

- QoS clamp measurements come from default, utility, and background taskpolicy launches.
- Core-class shares come from xctrace Time Profiler samples on Apple Silicon and are reported as sampled P-core vs E-core shares.
- The thread-scaling section is intentionally explicit that the current solver path is still single-threaded after the first CPU optimization pass.
