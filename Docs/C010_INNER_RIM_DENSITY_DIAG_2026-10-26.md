# C010 inner-rim density diagnosis at t = 6 (2026-10-26)

Authoritative fields: accepted t = 0 snapshot `/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core_refinement/fine_integrated/khout.000001.dat` and exact t = 6 snapshot `/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/varpress_w0_geom_cfl_long_20261026/C010/to6/khout.032479.dat`. Cell IDs were sorted and matched exactly. Radius bins use the t = 0 position around (12,12); every error below is against the same cell's initial value. Density and pressure changes are absolute fractional changes, mass weighted. Entropy is the same adiabatic invariant used by the existing C010 checkpoint audit.

| Initial radius | Cells | Mass-weighted abs density change | p99 abs density change | Mass-weighted abs pressure change | Mass-weighted abs entropy change |
|---|---:|---:|---:|---:|---:|
| 1–2 (inner taper) | 1,076 | 0.24684 | 1.29710 | 0.44701 | 1.44e-5 |
| 2–2.25 | 368 | 0.10488 | 0.45192 | 0.18072 | 1.60e-5 |
| 2.25–2.5 | 440 | 0.04010 | 0.18895 | 0.06774 | 2.71e-5 |
| 2.5–3 | 984 | 0.01119 | 0.05084 | 0.01864 | 5.85e-5 |
| 3–4 | 2,492 | 0.00199 | 0.00900 | 0.00321 | 1.93e-4 |
| 4–6 | 7,172 | 0.00027 | 0.00108 | 0.00030 | 1.68e-4 |
| 6–8 | 9,980 | 0.00004 | 0.00013 | 0.00004 | 2.77e-5 |

The current [density image](KEPLER_C010_DENSITY_T0_T6_2026-10-26.png) shows the same localization. The t = 6 whole-disk density error of 0.00621 hides a 0.10488 error in the first 0.25 radial units of the disk. Pressure changes there are close to the gamma = 5/3 adiabatic response to the density changes, while the entropy invariant barely changes. The independent ballistic-generator control below establishes the geometric origin of the disk's first-rim noise to high accuracy.

The operational priority remains to preserve C010 through its exact t = 10 gate and compare the next complete checkpoints. Do not infer long-time disk stability from this early t = 6 view (0.337 inner-edge orbits).

## Next complete checkpoint: t = 6.250024729493809

The native complete snapshot is `C010/to10/khout.036018.dat`, SHA-256 `4fdcbe0257bde7fd6cf7b8133eb786be24952fe5725be54c052cfba2775d1bc9`. The same accepted audit passes finite fields, exact individual masses, exact W = 0, positive internal energy, and positive volume. Global relative energy change is 5.3440e-10; minimum internal energy is 9.5934e-12.

| Initial radius | Mass-weighted abs density change at t=6.25 | Mass-weighted abs pressure change | Mass-weighted abs entropy change |
|---|---:|---:|---:|
| 1–2 (inner taper) | 0.26593 | 0.48468 | 1.51e-5 |
| 2–2.25 | 0.12554 | 0.21767 | 1.75e-5 |
| 2.25–2.5 | 0.05017 | 0.08518 | 2.73e-5 |
| 2.5–3 | 0.01320 | 0.02202 | 6.13e-5 |
| 3–4 | 0.00228 | 0.00369 | 2.02e-4 |
| 4–6 | 0.00029 | 0.00034 | 1.80e-4 |
| 6–8 | 0.00005 | 0.00005 | 3.03e-5 |

Whole-disk mass-weighted density error increased from 0.00620584 at t = 6 to 0.00747468 at t = 6.25, about 20.4%. The first 0.25 radial units of the disk increased from 0.10488 to 0.12554, about 19.7%. The small entropy change persists while density and pressure structure grows. Continue the chain and diagnose the geometric cause; no clamp or artificial damping follows from these measurements.

## Exact-orbit ballistic Voronoi control

The one-off [analysis script](ballistic_voronoi_control_t6_20261026.py) moved every initial generator to its exact circular-orbit position under the same softened potential at the native checkpoint time, without hydrodynamic force integration. SciPy then made an independent unweighted Voronoi diagram of all 71,080 generators. For the interior disk, its polygon areas agree with the native cell volumes at t = 6 and t = 6.25 to roughly 2.5e-8 mass-weighted relative error. The t = 0 interior areas match to the reported precision. The density prediction is the **same fixed cell mass divided by ballistic Voronoi area**; pressure prediction is its gamma = 5/3 adiabatic response.

| Time | Initial radius | Actual MW abs density change | Exact-orbit Voronoi prediction | MW abs residual | Signed-change correlation | MW abs pressure residual |
|---|---|---:|---:|---:|---:|---:|
| 6 | 2–2.25 | 0.1048818 | 0.1048950 | 3.82e-5 | 0.99999987 | 7.78e-5 |
| 6.2500247 | 2–2.25 | 0.1255348 | 0.1255535 | 5.33e-5 | 0.99999981 | 1.06e-4 |
| 6.5000792 | 2–2.25 | 0.1469215 | 0.1469476 | 7.68e-5 | 0.99999975 | 1.50e-4 |
| 6 | 2.25–2.5 | 0.0400977 | 0.0400994 | 9.56e-6 | 0.99999994 | 3.61e-5 |

Full machine-readable results: [t = 6](BALLISTIC_VORONOI_CONTROL_T6_2026-10-26.json), [t = 6.25](BALLISTIC_VORONOI_CONTROL_T625_2026-10-26.json), and [t = 6.5](BALLISTIC_VORONOI_CONTROL_T65_2026-10-26.json). This near equality is stronger than merely observing small trajectory errors: the observed inner-disk density and most pressure structure are quantitatively reproduced by the geometry of generators on perfect Kepler orbits. The native Voronoi tessellation agrees with the independent polygon-area calculation. The inner taper at r=1–2 has a weaker correlation (0.818 at t=6), so this conclusion is specifically strongest for the disk at r>=2.

For the disk, a Voronoi face moves differently from a material surface under differential rotation. Keeping each cell mass fixed while its geometric volume changes therefore creates density and adiabatic pressure changes even when all generators follow exact circular orbits. Intrinsic cell spin alone cannot alter an unweighted Voronoi cell's geometric volume. Under the user's no-steering and W=0 fixed-mass reduction constraints, the next distinct mechanism to test is **time-dependent Laguerre weights that control target cell volumes while retaining fixed mass**, with a consistent weighted-face velocity and conservation audit. Whether a fast stable weight evolution exists for this flow remains to be demonstrated; this diagnostic does not prove it.

## Accepted t = 6.500079194680843 checkpoint

The native complete snapshot `C010/to10/khout.039445.dat`, SHA-256
`821fb98951bb2e9afbd1448492996374b85ffd74625186d5a62261ce51b9076d`,
passes the finite, exact per-cell mass, W=0, positive internal-energy, and
positive-volume audit gates. Global relative energy change is `5.508e-10`,
and the smallest internal energy is `9.429e-12`. The whole-disk mass-weighted
density error grew from `0.0074747` at t=6.25 to `0.0089175` at t=6.5.
In the first disk ring, density change grew from `0.12553` to `0.14692` and
pressure change from `0.21767` to `0.25624`; its mass-weighted entropy change
is only `1.76e-5`. The independent exact-orbit Voronoi prediction remains
within `7.68e-5` mass-weighted absolute density residual there, with signed
per-cell correlation `0.99999975`. Thus the growing inner-rim signal is still
predominantly fixed-mass Voronoi geometry, not accumulated trajectory error.
