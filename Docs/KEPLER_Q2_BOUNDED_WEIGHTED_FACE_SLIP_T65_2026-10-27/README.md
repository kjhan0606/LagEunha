# Weight evolution and Kepler material-face speed at t=6.5

Independent CPU-only Slurm jobs **437356** (jitter) and **437357**
(aligned) completed `0:0` in 9 seconds each. The [auditor](../audit_kepler_bounded_weighted_face_slip_t65_20261027.py)
rebuilds the actual `t=6.5` ordinary Voronoi and solved target-volume
Laguerre diagrams. Across the same 145,671 faces shared by two initial
disk cells, it compares geometric face-normal speed with the **analytic
softened Kepler fluid velocity at that face centroid**. The weighted
cases use the accepted static weights and either zero, 16-CG-iteration,
or fully converged `Wdot` from the [rate probe](../KEPLER_Q2_BOUNDED_WEIGHT_RATE_T65_2026-10-27/README.md).

| Length-weighted mean absolute normal-speed difference | Aligned | 5% site jitter |
|---|---:|---:|
| Unweighted Voronoi, `W=0` | `9.66897e-5` | `1.24631e-4` |
| Static target-volume weights, `Wdot=0` | `9.68982e-5` | `1.25810e-4` |
| Target-volume weights, 16-iteration `Wdot` | `8.25383e-5` | `1.13839e-4` |
| Target-volume weights, converged `Wdot` | `8.22581e-5` | `1.13714e-4` |

The converged rate reduces this metric by about **14.9% aligned** and
**8.8% jittered** relative to W=0. Sixteen iterations are nearly identical
on this metric. Same-initial-shell faces remain important: their weighted
mean difference falls `1.11719e-4→1.00003e-4` aligned and
`1.13691e-4→1.00705e-4` jittered. Across different initial shells,
the corresponding values are `6.75299e-5→4.78299e-5` and
`1.45856e-4→1.38953e-4`.

The [aligned](aligned_slip.json) and [jitter](jitter_slip.json) receipts
report face counts, summed face lengths, RMS, and integrated absolute
speed differences for the whole disk and both face classes. They use the
same accepted snapshot and weight/rate hashes as the earlier geometry
audits; [file hashes](copied_inputs.sha256.json) verify the copies.

The volume-rate equation makes each cell's **signed sum** of geometric
face sweeps small; it does not set every face's normal speed equal to the
fluid speed. This measured residual means fixed mass plus volume control
alone is not a proof that Laguerre cells are material under Kepler shear.
The comparison uses the analytic benchmark velocity, not a native
face-reconstructed fluid state or actual conservative mass flux. Before
promoting a weighted hydro run, the operator must either establish a
material-face condition to sufficient accuracy or carry residual face
mass, momentum and energy exchange conservatively, while retaining the
`W=0` fixed-mass Voronoi reduction. These static results do not decide
the long-time stability or the necessary cell-spin contribution.
