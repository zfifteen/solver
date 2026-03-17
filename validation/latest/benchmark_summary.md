# Benchmark Summary

| Case | Reference Dataset | Primary Metric | Divergence L2 | Threshold | Status |
| --- | --- | --- | --- | --- | --- |
| couette_128 | analytic_couette_profile | profile_relative_l2_error=2.459100e-03 | 3.632610e-12 | profile_relative_l2_error <= 5e-3; divergence_l2 <= 1e-10 | pass |
| poiseuille_128 | analytic_poiseuille_profile | profile_relative_l2_error=2.239400e-05 | 3.930860e-14 | profile_relative_l2_error <= 5e-3; divergence_l2 <= 1e-10 | pass |
| lid_driven_cavity_re100_128 | Re100_centerline_literature_envelope | max_reference_relative_error=2.776500e-03 | 2.102650e-15 | named centerline sample points inside accepted envelope with miss <= 2%; divergence_l2 <= 1e-10 | pass |
| taylor_green_128 | analytic_taylor_green_decay_2d | normalized_energy_error=7.017980e-04 | 9.238060e-16 | normalized_energy_error <= 1e-2; divergence_l2 <= 1e-10 | pass |
| taylor_green_3d_64 | analytic_taylor_green_decay_3d | normalized_energy_error=2.488760e-03 | 8.727580e-16 | normalized_energy_error <= 1e-2; divergence_l2 <= 1e-10 | pass |
