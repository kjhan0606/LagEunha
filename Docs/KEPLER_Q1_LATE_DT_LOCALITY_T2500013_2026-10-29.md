# Late Q1 timestep locality and outer counterstreaming pairs

A read-only audit of the already saved native Q1 checkpoint at
t=2.5000127162090333 identifies where the very small stored timesteps occur.
No calculation was restarted, no new Slurm job was submitted, and the
ongoing halo40 trajectory and pending acceleration comparisons are unchanged.

The checkpoint is `diag_clip_swap_t2_to65_a100_20261028/khout.021137.dat`
under `/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1`.
Its SHA-256 is
`727508a1c46007126214ed73edf88e4725871ba76b7c2f3cab6b4e3c13938a26`.
The auditor checks the exact native record count/stride, complete unique IDs,
finite positive W0 state, unchanged masses and the matching step in the native
log. Stored per-cell dt is the last K4 force-stage value. It is not an RHS
recomputed on the accepted endpoint; accepted position and velocity fields
are used separately for the pair diagnostic below.

## What the checkpoint shows

| Quantity | Result |
| --- | ---: |
| Total cells | 147456 |
| Stored last-K4 dt minimum | 7.43675046e-5 |
| Stored last-K4 dt median | 0.2636861 |
| Cells below twice the minimum | 8 |
| Cells below four times the minimum | 35 |
| Initial-disk cells, 0.5<=r0<2 | 48680 |
| Minimum stored dt in that disk | 0.04764037 |
| Disk minimum / global minimum | 640.6073 |

The matching step 21137 accepted dt=0.00015575629367958754. The last four
printed CUDA minima before this record are 7.926e-5, 1.504e-4, 1.558e-4 and
7.437e-5. The native logger emits this record every fourth force call
(`Exam/exam_gpu_extract.c`), so these are K4 minima from four successive
steps, not the four stages of one step. The stored minimum agrees with the
last printed K4 value at its printed precision. That last-K4 minimum must
not be relabeled as the timestep already used by the accepted step.

Correction recorded on 2026-10-29: the earlier description called these
four records force-stage minima of one step. The numerical receipt and
the localization below are unchanged.

All eight cells below twice the minimum initially lie outside r0=2.75.
They form four mutually nearest pairs in the current periodic geometry:

| Generator IDs | Separation | Normal approaching speed |
| --- | ---: | ---: |
| 62986, 4979 | 7.5844131e-5 | 0.7843473 |
| 4059, 142835 | 7.5846734e-5 | 0.7843499 |
| 84469, 142476 | 7.5849198e-5 | 0.7843444 |
| 143396, 4620 | 7.5849268e-5 | 0.7843424 |

For every listed cell, the auditor searches all 147456 current generators
using minimum-image distances, rather than only testing the eight candidates.
With the input Courant coefficient 0.4, the acoustic estimate

    dt_pair = 2 C d / (cs_i + cs_j + max(0, closing_speed))

gives approximately the stored K4 dt: the estimate/stored ratio is
1.00000105–1.00000112. The estimate uses accepted fields and an EOS sound
speed, whereas the stored value came from the K4 predictor. This strong
agreement identifies a concrete candidate acoustic bottleneck, but does not
independently prove the winning native face CSR, every other limiter, or
the origin of every small dt during the subsequent evolution.

One example makes the counterstreaming visible without a mesh reconstruction.
Cell 4979 starts at (5.8046875,0.1953125), r0=3.9664271, and is at
(0.74358154,1.23010848), r=2.8677413, in the saved checkpoint. Its specific
angular momentum about (3,3) is now -0.2875029. Its nearest neighbor 62986
starts at r0=2.8683049 and retains specific angular momentum +1.6939713.
Their relative speed is about 0.95483, versus a summed sound speed about
0.03154. These are exterior counterstreaming pairs, not close pairs in the
initialized 0.5<=r0<2 disk.

This agrees with the previously identified Cartesian-periodic exterior
problem: circular trajectories outside r=3 encounter the box seams, while
the point-center gravity is not a smooth periodic continuation. The present
audit localizes one actual timestep bottleneck; it is not a new proof of
the entire crossing history or a claim that all later bottlenecks are exterior.
See also the existing periodic-seam and GIZMO comparison contracts.

## Interpretation for acceleration

Freezing this one stored dt distribution gives cell-update arithmetic ratios
of 55.06 for ideal independent updates and 37.19 for power-of-two bins.
These ratios exclude global mesh rebuilds, communication, additional timestep
bounds, changing limiters and conservative asynchronous face work. They are
not measured runtime gains and do not justify simply dropping global work.

The actionable evidence is that removing disk-internal shear constraints
alone would not remove this checkpoint's global bottleneck. Kernel/padding
acceleration, boundary/IC consistency and eventual conservative individual
stepping address different costs. Keep the current native trajectory and its
t>6 gates intact; assess any production change using actual pending comparison
results. A paper-equivalent long-time claim still requires the existing
GIZMO initial-condition and boundary contract, not only a passing active-disk
error norm in this different periodic exterior.

Artifacts: [reproducible reader](audit_kepler_late_dt_locality_20261029.py),
[full numerical receipt](KEPLER_Q1_LATE_DT_LOCALITY_T2500013_2026-10-29.json).
The bounded audit performs no tessellation or trajectory evolution. Native
cell spin, dynamic weights, t>6 and the GIZMO long endpoint remain unverified.
