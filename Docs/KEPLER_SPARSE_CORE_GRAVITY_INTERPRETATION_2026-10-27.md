# Central-force consistency of the corrected sparse Kepler IC

The [deterministic IC quadrature](audit_kepler_sparse_core_gravity_20261027.py)
examines the first ring of eight sites at `r=0.10080469` in each corrected
54,312/55,931-site sparse q2 glass. It constructs their exact initial
Voronoi polygons from the inner `r<0.4` generators, checks every polygon
vertex against excluded generators, and verifies the polygon areas against
the independently certified periodic sidecars. The maximum area disagreement
is `1.24e-12` relative. Product Gauss quadrature orders 32 and 64 agree in
the profile-weighted force ratio to better than `1.4e-12`. All inputs are
SHA-verified. The complete per-site [receipt](KEPLER_SPARSE_CORE_GRAVITY_IC_2026-10-27.json)
is preserved.

| Initial cells | Mean area-averaged radial force / generator force | Mean profile-weighted radial force / generator force | Mean cell-averaged profile density / generator profile density |
|---|---:|---:|---:|
| Sparse narrow, first 8 cells | 1.97757 | 1.30793 | 1.16757 |
| Sparse wide, first 8 cells | 1.97757 | 1.30793 | 1.16757 |
| Full q2, 8 nearest cells to the sparse sites (`r=0.110485`) | 1.00351 | 0.99333 | 1.00977 |

The force is the softened analytic acceleration
`g(x)=-x/(r²+0.01175²)^(3/2)`. The area average holds each finite-volume
cell's density constant; the profile-weighted average instead uses the
continuous target `rho(r)=0.01+(r/0.5)^3` within the core. The native RK4
driver applies `g` at each generator position. Its initial circular speed
balances that **point** acceleration. Thus the result quantifies a spatial
representation mismatch: the large first-ring cells do not resolve the
central force or density profile as well as the full q2 cells. The present
point-sampled assignment `m_i=rho(x_i)A_i` underestimates the target-profile
cell integral in these eight cells by about 16.8% relative to their assigned
mass. The eight assigned masses sum to `0.00131833`; the integrated target
profile sums to `0.00153925`, a difference of `0.00022092` or only
`1.86e-5` of the `11.904979` disk mass. This is a local resolution warning,
not a large whole-disk mass error. Self-gravity is disabled in this test.

This calculation does **not** establish an unstable native orbit, numerical
entropy production, or propagation into the disk. Area-averaged gravity by
itself would not be a self-consistent source replacement: the cell's
subcell velocity field and moving-face transport also contribute to the
finite-volume momentum balance. We therefore do not add a costly per-step
cell integration or change the submitted corrected IC. The result sets an
explicit core-resolution limitation for interpreting the sparse speedup.

Both corrected sparse native preflights subsequently passed the independently
certified periodic-area and fixed-mass gates at exact `t=.001`. An ID-matched
read of their first and last snapshots finds mean absolute generator radial
speed `1.40e-9` in the eight first-ring cells (`.1<=r<.2`); the full-q2
cells in that radial bin have `2.31e-10`. These tiny early generator speeds
show that the area-integrated force mismatch does not directly imply an
immediate orbit failure in the implemented point-source scheme. They do not
test the subcell momentum balance or later entropy production. At this
snapshot spacing, entropy changes are close to the stored single-precision
roundoff, so their relative sizes are not a reliable heating comparison.

Native comparisons at the same physical time must report the core and the
disk separately, including radial motion, pressure/entropy, and fixed-mass
cell-centre crossings of the inner disk edge. Those crossings are a transport
proxy, not an integrated fluid flux through the circle. The full q2 run is
the essential reference; a
passing sparse test alone cannot establish a general Kepler solver gain.
