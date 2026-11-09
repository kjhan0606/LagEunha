# Original nine 1D tests with shared material-face PLM

At the user's request, stop deep investigation of the isolated Blast peak and execute the full original 1D suite with the current candidate. All nine now reach their original target times under one fixed-mass PLM/minmod/SSPRK2 operator. No case-specific viscosity, entropy fix, floor or energy reset. Original initial states/gamma/physical intervals/endpoints and wall boundaries are taken from the retained fv_original_cases.h. Sod uses its retained reflecting boundaries. Shu uses 800 cells in the physical interval; all other cases use200. Open problems retain signal*time guard padding and moving transmitting outer boundaries. Total cell counts are reported separately.

This suite uses original moving-contact speed1 (the preliminary four-case comparison used .5). Initial conserved amounts are integrated on each real material interval, including the analytic sin(5*pi*x) density integral in Shu. Guard independence has not been separately established for this new material operator.

| Case | Physical N | Target | Density average L1 |
|---|---:|---:|---:|
| sod | 200 | 0.2 | 0.6883% |
| blast | 200 | 0.038 | 4.998% |
| noh | 200 | 0.6 | 0.4029% |
| shu | 800 | 1.8 | 0.5933% |
| lax | 200 | 0.13 | 0.563% |
| doublefan | 200 | 0.15 | 2.546% |
| collision | 200 | 0.035 | 0.4105% |
| contact_open | 200 | 0.25 | 2.59e-11% |
| moving_contact_open | 200 | 0.25 | 1.122e-11% |

All original cell masses are bit-identical, interfaces exactly shared, all final cell widths and pressures positive. Whole computational-domain energy balance includes boundary pressure work; largest normalized error <4e-15. Wall boundary impulse is included in momentum balance. These conservation checks do not imply pointwise accuracy.

Known residuals: Blast maximum density8.152 vs independent reference6.442; Noh wall-adjacent maximum4.494 vs finite-pressure exact3.99999; collision maximum31.494 vs31.043. Double rarefaction retains an under-resolved central region (velocity relative L1 approximately14.93% within[0,1], pressure5.38%). Its density error2.546% in[0,1] is not directly comparable with the earlier1.019% over the larger moving computational interval. No failed result or residual peak is hidden.

References: independent Athena analytic Riemann utility for Sod/Lax/collision/Noh, analytic symmetric rarefaction, analytic constant-pressure contacts; retained independent AthenaN6400 Blast andN12800 Shu at full target time. Numerical references are not exact. Accuracy integration uses actual material-cell intervals clipped to original physical domains. Plots show actual horizontal cell-average segments plus centre markers, not smoothing.

Density, pressure, velocity PNG/PDF figures, snapshots, case specifications, run receipts and source hashes are stored here. Native C simulations use one CPU per process, three independent cases concurrently. All nine simulations have completed; no background computation is claimed. Existing production OpenMP/MPI/CUDA code and previous outputs unchanged. This is1D evidence, not multidimensional convex Laguerre/rotation or long-Kepler validation.

Evidence directory: /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_full_original_nine_20261105_v1
