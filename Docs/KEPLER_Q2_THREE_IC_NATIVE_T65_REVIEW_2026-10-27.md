# Three bounded-q2 Kepler ICs at native t=6.5

The [hash-verified three-IC comparison](KEPLER_Q2_THREE_IC_NATIVE_T65_2026-10-27.json)
combines the accepted aligned, ring-coherent 5%, and independent site-angle
jitter 5% W=0 CUDA runs at exact `t=6.5`. All three endpoint auditors pass
17/17 checks, including fixed masses, positive states, exact bounded-area
tiling, no injected energy, and orbital/energy gates. The ring-coherent job
**437389** completed `0:0`; its audited endpoint SHA-256 is
`d911d4e122aed677b2431fa7e7a3890c87d9a372b06c1d6abd364188bc0fb554`.
The [common-scale native density figure](KEPLER_Q2_THREE_IC_NATIVE_DENSITY_T65_2026-10-27.png)
shows the surface-density fields above and relative density changes below;
nearest-site coloring is for display, while all tabulated metrics use the
original cells and masses.

| Mass-weighted disk measure at t=6.5 | Aligned | Ring-coherent 5% | Independent site jitter 5% |
|---|---:|---:|---:|
| Mean absolute density change | 1.86493e-4 | 2.06530e-4 | 1.69828e-3 |
| Exact-orbit Voronoi density change | 1.85727e-4 | 2.05525e-4 | 1.47879e-3 |
| Native-minus-exact-orbit density residual | 2.88254e-6 | 3.75748e-6 | 3.44253e-4 |
| Mean circular-orbit position error | 1.69938e-7 | 1.86478e-7 | 5.21458e-6 |
| Mean absolute entropy-log change | 7.47903e-6 | 8.18232e-6 | 1.49389e-5 |

The native ring-coherent density measure is only **1.11×** aligned, while
independent site jitter is **8.22×** ring-coherent. The native-minus-ballistic
residual is **91.6×** greater for site jitter than ring-coherent, and its
orbit error is **28.0×** greater. This directly supports the hypothesis that
loss of regular spacing **within each ring**, rather than a small relative
angular phase between neighboring rings, drives most of this IC-sensitive
error at the tested amplitude and resolution. The ring-coherent exact-orbit
geometry prediction was also quantitatively accurate for its native density
change. The separate [signed entropy audit](KEPLER_Q2_THREE_IC_SIGNED_ENTROPY_T65_2026-10-27.json)
finds mass-weighted mean changes of `−4.94e-6` aligned, `−5.15e-6`
ring-coherent and `−3.37e-7` site-jitter, versus the larger absolute local
changes shown in the table. The site-jitter result does not show a net
positive entropy drift at this time. Signed cancellation cannot exclude
local dissipation, but the observed density noise alone should not be called
irreversible heating.

The [initial pressure-force series](KEPLER_Q2_INITIAL_PRESSURE_ONSET_REVIEW_2026-10-27.md)
adds a time ordering: all three ICs have zero pressure acceleration at
`t=0`, but by `t≈0.825` independent jitter has `15.5×` the aligned
pressure-acceleration RMS. Its density change at that time is already
almost entirely predicted by exact-orbit Voronoi geometry. This supports
geometry-generated pressure contrast as an early feedback trigger, while
the integrated cause of the eventual orbit error remains to be isolated.
The separate [same-t0 pressure-force-off control](KEPLER_Q2_EARLY_PRESSURE_CAUSAL_TRIAD_T0825_2026-10-27.md)
finds that almost all **t=0.825** orbit deviation vanishes without the
pressure force, even though the geometric density pattern remains.

The [exact-orbit W=0 geometry reference through t=600](KEPLER_Q2_THREE_IC_BALLISTIC_TO600_2026-10-27.json)
keeps the same initial masses and moves only generators on analytic softened
circular orbits. Its disk mean absolute density change remains near
`(1.9–2.3)e-4` aligned, `(2.1–2.3)e-4` ring-coherent and
`(1.47–1.49)e-3` independent site-jitter from `t=6.5` through `t=600`.
Thus exact orbital shear of the irregular Voronoi geometry creates a
persistent noise level but does not by itself cause runaway growth in this
read-only construction. The existing independent-jitter native `t=10`
endpoint reaches `2.37446e-3`, with native-minus-exact-orbit residual
`1.19668e-3`; the extra growth requires a dynamical response beyond the
kinematic area change. The exact-orbit construction is a reference, not a
strict lower bound for every possible numerical trajectory and not a native
`t=600` stability result.

Each IC uses its own Voronoi area to set initial fixed cell mass from the same
continuum surface-density profile. The cells therefore do **not** have
identical discrete masses across the three runs. This is a legitimate
mechanism comparison, but polar order, sparse-core force resolution and the
fixed polygonal wall remain special. A stable aligned or ring-coherent disk
cannot establish general Kepler accuracy or equivalence to the long GIZMO
benchmark. Follow with a nonaligned resolution/timestep study and boundary
sensitivity, using density, orbit, native-minus-ballistic and signed entropy
together. Increasing the weighted Jacobi iteration count from 16 to 64 did
not remove the independent-jitter `t=10` density failure, so more linear
iterations alone are not the next physical remedy.
