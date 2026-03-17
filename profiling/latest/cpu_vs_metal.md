# CPU vs Metal Taylor-Green Comparison

These comparisons use the same Taylor-Green executable and configuration surface, with the Metal path selected via `--backend metal`.

The current Metal slice is velocity-first: it uses float working storage on device and a final CPU projection cleanup before result publication, so velocity / energy agreement is the primary correctness comparison for this milestone.

The `256^3` smoke case is intentionally metrics-only here to avoid materializing giant ASCII VTK snapshots during a long profiling run.

| Case | Comparison Mode | CPU Status | Metal Status | CPU Seconds | Metal Seconds | Speedup | Velocity Rel L2 | Pressure Rel L2 | Energy Rel Diff | CPU Div L2 | Metal Div L2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| taylor_green_3d_64 | full_vtk_diff | pass | pass | 2.680 | 1.930 | 1.389x | 1.473779e-03 | 1.412473e-02 | 7.239673e-05 | 8.732090e-16 | 3.269290e-16 |
| taylor_green_3d_256_smoke | metrics_only | skipped | skipped | 33.160 | 15.020 | 2.208x | n/a | n/a | 0.000000e+00 | 1.305000e-15 | 1.305610e-15 |
