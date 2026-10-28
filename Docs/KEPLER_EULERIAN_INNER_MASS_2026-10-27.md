# Kepler inner-circle mass from native Voronoi cells (2026-10-27)

The existing `r=0.5` crossing count follows **cell generators**, not the gas
mass occupying the fixed Eulerian circle. The
[offline auditor](audit_kepler_eulerian_inner_mass_20261027.py) reconstructs
each W=0 Voronoi polygon, intersects it exactly with the circle, and sums
`m_i * intersection_area / polygon_area`. This is a piecewise-constant cell
reconstruction that conserves each complete cell's fixed mass. It does not
integrate face flux over time or reconstruct subcell density gradients.

The candidate search is certified: for any point within radius `R`, a core
site at radius `r_min` is within `R+r_min`; a site farther than
`2R+r_min` from the center cannot own that point. Every candidate polygon
must be bounded, have positive area, agree with native volume to `1e-4`
relative, and satisfy `edge_distance > 2*circumradius` so omitted periodic
images cannot change it. The circle-polygon formula passed analytic square,
half-disk, and enclosed-square checks to `1e-12` and 25 random convex
polygon comparisons against a 4,096-segment circle to `1.93e-8` absolute.

Independent grammar CPU Slurm job **1223431** completed all three native
cases in 24 s, exit `0:0`, peak batch RSS 540 MB. The staged inputs and
script hashes are at
`/gpfs/kjhan/LagForce/Kepler_InnerMass_Preflight_20261027/inputs.sha256`.
Individual receipts are [full q2 at t=0.2062646](KEPLER_EULERIAN_INNER_MASS_FULL_Q2_T020626_2026-10-27.json),
[corrected sparse narrow at t=.001](KEPLER_EULERIAN_INNER_MASS_SPARSE_NARROW_T0001_2026-10-27.json),
and [corrected sparse wide at t=.001](KEPLER_EULERIAN_INNER_MASS_SPARSE_WIDE_T0001_2026-10-27.json).

| IC | Initial `M(r<0.5)` | Difference from analytic target | Relative target error | Native change at audited time | Generator-crossing change |
|---|---:|---:|---:|---:|---:|
| Full q2 | 0.321615283 | −0.000397964 | −0.124% | −9.42e-9 at t=.2062646 | 0 |
| Sparse narrow | 0.324132920 | +0.002119673 | +0.658% | −1.90e-11 at t=.001 | 0 |
| Sparse wide | 0.324086038 | +0.002072791 | +0.644% | −1.91e-11 at t=.001 | 0 |

The exact continuum target from native `rho=0.01+(r/0.5)^3` for `r<0.5`
is `pi*(0.5)^2*(0.01+2/5) = 0.322013246993`. Although all three ICs have
identical generator positions in the main `0.5<=r<2` disk and the same disk
mass, sparse narrow/wide begin with **0.783%/0.768% more mass inside the
inner circle than full q2**. The sparse mass assignment is
`rho(generator)*Voronoi_area`; the coarser core and transition cells sample
the target profile differently. This is a pre-existing IC difference, not
mass created during the short native runs.

The geometric boundary correction `Eulerian M(<.5) - mass of generators
inside .5` is `0.008940157131` for all three ICs to `5e-13`. Consequently,
the sparse/full enclosed-mass differences agree with their core-generator
mass differences to `5e-13`: the measured initial bias comes from the changed
`r<.5` mass sampling, rather than a different tessellation at the disk edge.
This exact cancellation is consistent with the 49,038 bitwise-identical disk
generator positions established in the whole-disk geometry control.

At these audited times, native inner mass barely changes, so the zero
generator-crossing proxy happens to agree with the negligible evolution.
Its absolute initial enclosed-mass estimate is low by about `0.00894`
because many cells straddle the circle; zero future crossings therefore
cannot by itself certify zero Eulerian transfer. The largest candidate
native/Voronoi area mismatch was `5.89e-8`, and the minimum omitted-image
safety ratio was `14.38` (sparse) or `81.51` (full endpoint).

For `t=.825`, `t=10`, and later native checkpoints, repeat this audit on
each accepted state. Compare **changes from each IC's own initial enclosed
mass**, and report the initial bias separately. Keep the submitted sparse
continuations as controls. A new IC that assigns target-profile cell mass
by one-time polygon quadrature could reduce this bias while retaining fixed
mass and the same GPU evolution cost; it needs a separate preflight and
same-time native comparison before replacing the current IC.

## Native q1 halo-12 endpoint at exact t=0.825

The saved q1 recovery checkpoint `khout.010856.dat` completed native
evolution at exact `t=.825` with no clipped halo cells, nonpositive internal
energy, floor injection, or dual-energy injection. Its Slurm batch returned
`FAILED 1:0` only after the native solver exited zero: the staged radial
auditor lacked its Python import. The state itself passes the independent
[strict 15-check gate](KEPLER_Q1_HALO12_T0825_STRICT_2026-10-27.json)
and [radial audit](KEPLER_Q1_HALO12_T0825_RADIAL_2026-10-27.json)
in grammar CPU job **1223461**; no GPU calculation was repeated.

The [boundary-aware circle integral](KEPLER_EULERIAN_INNER_MASS_Q1_T0825_2026-10-27.json)
in CPU job **1223468** used `np.hypot`, matching the radial auditor's
generator-radius convention. It finds:

| Diagnostic | t=0 | t=.825 | Change |
|---|---:|---:|---:|
| Eulerian piecewise-constant mass inside r=.5 | 0.322388456934 | 0.322388467914 | +1.09805e-8 |
| Mass of generators at r<.5 | 0.320452050518 | 0.324325784984 | +0.003873734466 |

The q1 IC places **203 generators**, carrying `0.04955853` mass, within
`1e-12` of the `r=.5` circle. A strict center-inside test is therefore
ill-conditioned: even using `np.linalg.norm` instead of `np.hypot` moves
three initial sites, mass `0.0007333672`, across its Boolean threshold.
The fixed-circle mass integral has no such point-classification jump; all
contributing Voronoi polygons are certified independent of omitted periodic
images (minimum safety ratio `91.27` at the endpoint), and native/Voronoi
areas agree to `5.72e-8` relative. The apparent `0.00387` inward cell-center
mass is **not** evidence of comparable physical inflow. It entirely reflects
crossings of generators initialized on the diagnostic boundary. Future
transport conclusions must use the Eulerian integral and, where needed,
face-based time-integrated flux rather than this center-crossing proxy.

An ID-matched check of the two native snapshots confirms that **all 108**
reported generator crossings (62 inward, 46 outward; total crossed mass
`0.02636332`) originate from those `r=.5` initial sites. No generator
starting farther than `2.78e-16` from the boundary crosses it, and the
largest radial displacement among crossing sites is only `4.11e-8`. This
directly identifies the apparent transport as a boundary-classification
artifact in q1.
