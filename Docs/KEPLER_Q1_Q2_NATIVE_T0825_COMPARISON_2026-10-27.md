# Matched native q1/q2 polar Kepler comparison at t=0.825

The accepted CUDA+MPI+OpenMP W=0 fixed-mass q1 and q2 runs both have a
complete exact `t=.825` snapshot. The q1 [strict gate](KEPLER_Q1_HALO12_T0825_STRICT_2026-10-27.json)
and new q2 [strict gate](KEPLER_Q2_HALO12_T0825_STRICT_2026-10-27.json)
each pass all 15 checks: positive state, exact mass and W=0, periodic area
sum, `rho=m/V`, zero floor/dual-energy injection, orbital and density
thresholds, disk angular momentum, and whole-box energy. The q2 native
solver and Slurm job **436995** completed with exit `0:0`; its final radial
audit also passes. The [machine comparison](KEPLER_Q1_Q2_NATIVE_STRICT_COMPARISON_T0825_2026-10-27.json)
requires both passing receipts at identical physical time.

| Initial `0.5 <= r < 2` disk | q1 | q2 | q2/q1 |
|---|---:|---:|---:|
| Native mean absolute relative density change | 5.23520e-4 | 1.16521e-4 | 0.2226 |
| Exact-orbit Voronoi density change | 5.23767e-4 | 1.16545e-4 | 0.2225 |
| Native/ballistic pointwise density residual | 4.80819e-7 | 8.76967e-8 | 0.1824 |
| Absolute entropy-log change | 1.53318e-4 | 1.29736e-4 | 0.8462 |
| Signed entropy-log change | -8.45393e-6 | -1.89861e-5 | — |
| Exact-orbit position error | 1.11215e-8 | 2.50253e-9 | 0.2250 |

The independent q2 [native-versus-orbit geometry receipt](KEPLER_Q2_NATIVE_BALLISTIC_T0825_2026-10-27.json)
comes from grammar CPU Slurm job **1223641**, completed `0:0` in 23 seconds.
Its inputs and reader layout are hashed in
`/gpfs/kjhan/LagForce/Kepler_Q2_BallisticResidual_T0825_20261027/inputs.sha256`.
The q2 pointwise native/ballistic density residual is only **0.0753% of
the native absolute density-change norm**. This is a comparison of two
positive norms, not an exact signed error decomposition. The corresponding
q1 fraction was 0.0918%. The native and exact-orbit disk density patterns
match visually on the same color scales in the [q1 density field](KEPLER_POLAR_Q1_NATIVE_DENSITY_T0825_2026-10-27.png)
and [q2 density field](KEPLER_POLAR_Q2_NATIVE_DENSITY_T0825_2026-10-27.png).
The top panels show actual native density; the lower panels show absolute
relative density change and the exact-orbit Voronoi prediction.

The density gain is much larger than the entropy gain. In the initial
`0.75 <= r < 1` band, q2's absolute entropy-log error is **1.403 times
q1's**, even while its density error is only 0.238 times q1's. The signed
disk-wide entropy changes are negative in both layouts; neither early
snapshot demonstrates net irreversible heating. A future cell-spin or
Laguerre modification must be judged against native entropy/torque and
orbital evolution rather than against the density image alone.

Both layouts use the same outer Cartesian sites. Their initial `r>=2.75`
cells have **identical** mean absolute density change `0.0922209` and mean
absolute radial speed `0.0715428` at `.825`, from the q1 and q2
[radial receipts](KEPLER_Q1_HALO12_T0825_RADIAL_2026-10-27.json)
and [q2 radial receipt](KEPLER_Q2_HALO12_T0825_RADIAL_2026-10-27.json).
Thus a quieter disk mesh has not repaired the non-equilibrium periodic
outer boundary. The q1/q2 disk site distributions also differ, so this is
an IC comparison, not a one-operator ablation.

At `t=.825`, the q2 polar IC preserves circular orbits and density more
accurately than q1 in native gas evolution; the early density error is
mostly predicted by fixed-mass Voronoi geometry. The result does **not**
prove suppression of late shear heating or survival past the previous
`t=6–10` failure region. The queued native q1/q2 to `t=6.5` and corrected
sparse to `t=10` remain the decisive next physical gates. A paper-level
comparison still requires the native density field at `t=120` and the
published longest plotted time `t=600`.
