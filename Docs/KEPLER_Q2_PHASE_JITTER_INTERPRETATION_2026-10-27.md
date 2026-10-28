# Kepler q2 polar IC: phase-jitter geometry control (2026-10-27)

The [Astra IC review](KEPLER_POLAR_NATIVE_CONTROL_2026-10-26.md#ic-admissibility-and-benchmark-interpretation-astra-review-2026-10-27)
identified a missing test: does q2 work only because its sites are precisely
aligned on each circular ring? The [read-only control script](audit_kepler_q2_phase_jitter_geometry_20261027.py)
starts from the complete native q2 t=0 generators (147,456 sites, 49,038
disk sites). In `0.48<=r<2.02`, it offsets each site angle by a fixed random
fraction of its ring's `2π/N` spacing, preserving its radius. Two controls
use the same amplitude: independent site offsets, or a common offset for
every site on a given ring. The disk's initial mass is proportional to its
*perturbed* initial Voronoi area, so `rho=m/A` starts at the same 1.01
plateau. Every site then follows its exact softened circular orbit; the
Voronoi areas are recomputed at t=6.5 and t=600. All tested disk polygons
pass the omitted-periodic-image safety test by a factor greater than 40.

The table gives disk area-weighted mean `|A(0)/A(t)-1|`, equal to fixed-mass
relative density drift for the 1.01 plateau. One fixed random seed
(`20261027`) was used for each perturbation mode.

| Perturbation | Fraction of ring spacing | t=6.5 | t=600 |
|---|---:|---:|---:|
| None | 0 | 0.00018573 | 0.00019699 |
| Each ring rotates coherently | 0.05 | 0.00020553 | 0.00021135 |
| Each ring rotates coherently | 0.10 | 0.00022680 | 0.00022254 |
| Each site has its own offset | 0.05 | 0.00148275 | 0.00146630 |
| Each site has its own offset | 0.10 | 0.00295217 | 0.00291816 |

Independent 5% site jitter raises the t=600 geometry error by about 7.4
times; coherent 5% ring phase jitter raises it by only about 7%. The
site-jitter error is still about 0.15% in the mean. The contrast shows that
along-ring regularity is important to q2's unusually small geometry error.
Radial-versus-tangential spacing alone is an incomplete explanation.

This is not a native gas calculation. The result does not determine entropy
production, face work, torque, positivity, or survival to t>6. One seed is
enough to expose sensitivity, not to quantify its distribution. A subsequent
native comparison should use the same q2 physics and solver with a small
site-jitter IC whose initial masses are recomputed from the native initial
tessellation; compare it with unperturbed q2 at equal physical time through
at least t=10. Keep the Cartesian control for the benchmark's generality.

Machine-readable results: [independent-site jitter](KEPLER_Q2_PHASE_JITTER_GEOMETRY_2026-10-27.json),
[coherent-ring jitter](KEPLER_Q2_RING_PHASE_JITTER_GEOMETRY_2026-10-27.json).

## Native matched control queued (2026-10-27)

The [5% independent-site jitter generator](make_kepler_q2_site_jitter_glass_20261027.py)
now transforms the byte-verified accepted full-q2 glass into
[a native glass](KEPLER_Q2_SITE_JITTER_005_2026-10-27.bin), SHA-256
`754c6c60acb516404f6799111dac6605ad58c04baa166f706144c046e6a5d080`.
Its [manifest](KEPLER_Q2_SITE_JITTER_005_2026-10-27.json) confirms the
same 147,456 sites and 49,038 disk sites, 50,481 perturbed disk/guard sites
on 70 rings, maximum radius roundoff `8.88e-16`, and minimum pair spacing
`0.00979103`. Nonpolar sites are byte-identical. The generator uses the
same global random sequence and angular rule as the offline geometric audit.

Native H200 preflight **437143** uses the accepted q2 CUDA+MPI+OpenMP binary
SHA-256 `7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`,
the same physical settings, and the same one-time `m=rho_target*V(0)`
initialization. Its only IC change is the glass. The [preflight runner](kepler_q2_site_jitter_preflight_20261027.slurm)
and actual staged inputs passed shell syntax and SHA checks; it is pending
GPU quota at reduced scheduling priority. It must pass exact `t=.001`
fixed-mass, W=0, positivity, initial-density and orbit checks.

Conditional H200 job **437145** has verified `afterok:437143`. Its
[runner](kepler_q2_site_jitter_to020625_20261027.slurm) and
[checkpoint preparer](prepare_kepler_q2_jitter_to020625_20261027.py)
refuse an unaudited or incomplete preflight and restart the same binary.
The endpoint is the **actual** accepted q2 saved time
`t=0.20626460064652694` (`khout.002700.dat`), not its nominal `.20625`
dump clock. The [native matched auditor](audit_kepler_q2_jitter_t020625_20261027.py)
will compare fixed mass, W=0, positive state, signed/absolute entropy and
density drift, radial velocity, angular momentum, and ballistic geometry
residuals at this identical physical time. Both jobs use H200:1, 32 GB and
positive Slurm nice so they do not outrank the existing baseline/sparse
sequence. **No native jitter outcome exists yet.**
