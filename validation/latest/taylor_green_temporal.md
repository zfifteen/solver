# Taylor-Green Temporal Refinement

Temporal order is measured with a self-convergent final-field comparison on the same 128 x 128 grid.

| CFL | dt | Normalized KE Error | Velocity Relative L2 | Divergence L2 | Self Difference | Order |
| --- | --- | --- | --- | --- | --- | --- |
| 1.00 | 2.272730e-02 | 7.018390e-04 | 2.591940e-03 | 3.642860e-15 | 3.406626e-05 |  |
| 0.50 | 1.190480e-02 | 7.017960e-04 | 2.587030e-03 | 1.409960e-15 | 8.858612e-06 | 1.943 |
| 0.25 | 6.097560e-03 | 7.017980e-04 | 2.585750e-03 | 9.238060e-16 |  |  |
