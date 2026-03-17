# CPU vs Metal Taylor-Green Comparison

These comparisons use the same Taylor-Green executable and configuration surface, with the Metal path selected via `--backend metal`.

The current Metal slice is velocity-first: it uses float working storage on device and a final CPU projection cleanup before result publication, so velocity / energy agreement is the primary correctness comparison for this milestone.

| Case | CPU Status | Metal Status | CPU Seconds | Metal Seconds | Speedup | Velocity Rel L2 | Pressure Rel L2 | Energy Rel Diff | CPU Div L2 | Metal Div L2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| taylor_green_3d_64 | pass | pass | 2.990 | 2.550 | 1.173x | 1.473778e-03 | 1.412472e-02 | 7.239673e-05 | 8.732090e-16 | 3.290580e-16 |
| taylor_green_3d_256_smoke | skipped | skipped | 49.760 | 27.710 | 1.796x | 1.841304e-05 | 3.797453e-03 | 0.000000e+00 | 1.305000e-15 | 1.304940e-15 |
