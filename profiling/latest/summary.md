# Performance Summary

## Hardware

- CPU: Apple M1 Max
- Performance cores: 8
- Efficiency cores: 2
- Memory: 32.0 GiB

## Findings

- Recommended default execution mode: benchmark build profile with the default unclamped scheduler policy. It was the fastest measured policy on this machine.
- Current compute-thread recommendation: 1. The solver path remains single-threaded after the first M11 optimization pass, so thread scaling is still a baseline-only study for now.
- Advection microbenchmark throughput: 1.328590e+07 cell updates/s with a 0.852 GB/s lower-bound bandwidth estimate.
- Pressure microbenchmark throughput: 1.668890e+06 unknown updates/s with 16.00 average iterations.
- End-to-end Taylor-Green benchmark throughput: 1.033452e+06 cells/s.
- Targeted Metal backend result: 1.173x faster than benchmark-profile CPU on `taylor_green_3d_64`, with velocity relative L2 drift 1.473778e-03.
- Targeted Metal backend smoke throughput: 1.796x faster than benchmark-profile CPU on `taylor_green_3d_256_smoke`.
- Metal correctness in this first slice is tracked primarily through final velocity / energy agreement plus the analytic benchmark gate, because the device path uses float working storage before a final CPU projection cleanup.
- Default-policy hotspot categories were dominated by pressure-solve, predictor/ADI, and advection work, with the top sampled category counts: {'other': 40, 'pressure_solve': 25, 'predictor_adi': 17}.

## Comparison To Milestone 10 Baseline

- Advection microbenchmark throughput: 1.046x vs the stored Milestone 10 M1 Max baseline.
- Pressure microbenchmark throughput: 5.655x vs the stored Milestone 10 M1 Max baseline.
- Taylor-Green end-to-end throughput: 3.554x vs the stored Milestone 10 M1 Max baseline.

## Notes

- QoS clamp measurements come from default, utility, and background taskpolicy launches.
- Core-class shares come from xctrace Time Profiler samples on Apple Silicon and are reported as sampled P-core vs E-core shares.
- The thread-scaling section is intentionally explicit that the current solver path is still single-threaded after the first CPU optimization pass.
