# Q1 t=3 radial failure: exterior generators compress disk cells

## Current execution state

Native CUDA+MPI+OpenMP job **438957** reached exact t=3, step 27642 and exited its solver with code 0. The 15 native checks pass, but the additional radial continuation gate passes only 8/9: initial-radius 1<=r0<2 density change is **0.0342085884**, exceeding the unchanged 0.01 threshold. The batch therefore ended FAILED 1:0 and submitted no t=4 successor. The complete checkpoint remains preserved at:

`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/halo40_from_t2654918_a100_20261028/to3p0/khout.027642.dat`

SHA-256: `cfa09e4f0f8b074bed39b35fe15d2e48fb14d1f3bb7c86a3968aaf2539c02bae`.

Global energy drift from t=0 is -2.80713e-8; disk angular-momentum drift is 3.04139e-10. Floor, dual-energy and entropy-floor injection ledgers are zero. Disk mean absolute entropy-log change is 3.59648e-4, with signed mean -4.15756e-5; this is not evidence for net irreversible shear heating.

The existing seven-case post-t3 dispatcher was not reached because its prerequisite radial gate failed. The prepared halo40 promotion workflow remains uninstalled and unselected. Its old conditional approval does not apply to this failed campaign gate. No gate was relaxed and no physical trajectory was restarted.

## Independent saved-state diagnostic

CPU Slurm job **1234685** on grammar004 completed 0:0 in 54 seconds, using under 600 MiB peak RSS. It reads the original 384² initial state and the accepted t=2.5000127162, 2.6549183919 and 3 states. It evolves no fluid and writes one JSON receipt and one density PNG. Job 1234683 is retained as a superseded diagnostic: using `np.linalg.norm` instead of the production `np.hypot` classified two sites on exact radial boundaries differently. The final diagnostic matches all 48680 production disk sites and the production density norms exactly.

For actual generator coordinates, an independent SciPy Voronoi reconstruction reproduces the stored disk volumes within maximum relative **5.624e-8** at t=3. Every generator is checked inside the primary 6x6 box, and every selected disk polygon satisfies edge-distance > twice its circumradius (minimum ratio 45.79). This excludes every omitted periodic image from cutting these cells, so the reconstructed disk polygons are the full periodic result. Summed area alone was not used to claim individual volume accuracy.

Define exterior provenance by initial radius r0>=2.75. Direct support means an exterior generator actually shares a Voronoi face with a disk generator.

| Saved time | Disk cells with exterior support | 1<=r0<2 native density-change norm | Same inner positions, isolated exact-orbit exterior |
| --- | ---: | ---: | ---: |
| 2.5000127162 | 0 | 4.33602e-4 | 4.33602e-4 |
| 2.6549183919 | 384 | 3.06387e-3 | 4.07960e-4 |
| 3 | 5018 | 3.42086e-2 | 3.95918e-4 |

At t=3, directly affected disk cells contribute **98.312%** of the disk absolute density-change norm. Their exterior supporting sites originated at r0>=3.16237, outside the inscribed radius-3 circle of the periodic box. There are 1856 exterior-provenance generators now inside r<2. The largest disk density changes approach +99% even while the affected disk generator's orbit error can be about 1e-8.

Two exterior counterfactuals were computed while preserving every actual inner (r0<2.75) position: restore exterior coordinates to t=0, or place them on their isolated exact circular orbits. Both reduce the 1<=r0<2 density-change norm to 3.95918e-4. Replacing all positions with isolated exact orbits gives 3.98286e-4. The actual full-disk norm is 2.75373e-2; the exterior replacement gives 4.97811e-4. Thus this late increase is primarily caused by exterior support intrusion, rather than a failure to compute the actual native Voronoi areas or large disk orbit departure.

These replacements are diagnostic configurations, not evolved fluid states or alternate solutions to the same periodic problem. Their masses and energies must not be promoted into production. The isolated reference intentionally has different exterior boundary/support data. The signed density m=4 amplitude is 0.01897, consistent with the four-sector pattern, but Fourier symmetry alone is not the causal test; changing only exterior positions is the stronger evidence.

## Volume-independent cross-check

Sharp radial mass cuts coincide with exact polar rings, so tiny orbital errors can move many sites across an exactly aligned boundary. The supplementary reader uses tanh radial windows of width 0.02 and sums fixed masses, independently of cell volumes. For the smooth 0.5–2 window, mass changes by 0.00463533: 0.00463331 comes from exterior-provenance sites and only 2.02085e-6 from the remaining sites. This is not a face-integrated flux or a direct conversion of the cell-density norm into bulk compression.

## Consequences for the next development step

The q1 periodic chain has now demonstrated exterior contamination reaching the physical disk before t=6. Increasing halo size repairs search completeness but cannot repair this global boundary/initial-condition inconsistency. The current evidence supports prioritizing a consistent nonperiodic boundary in the user's native GPU solver, with explicit force/work and boundary impulse/energy/angular-momentum accounting. Spin, weights and individual stepping do not remove this identified cause. Preserve W0/fixed mass/no steering/Re=infinity for the boundary comparison.

The earlier closed polygon q2 control reached t=30 but eventually stopped at t=48.3753 on an exterior generator crossing. Its absence of a validated collision/impulse rule and the Rayleigh-unstable tapered-potential alternative are already documented; neither should be replayed as an assumed cure. Astra recommends completing a conservative slip-wall path in the existing bounded native solver: boundary-matched reflection, accepted-step event handling with rollback, and separate pressure/collision impulse and torque ledgers. The first native host/device impulse primitive is now implemented and tested; RK4 event integration remains outstanding. Use polygon normals for the existing polygon geometry; radial reflection is reserved for a truly circular tessellation boundary. The clean Q1 boundary IC must preserve the original disk, and boundary-radius/shape dependence must be checked before any GIZMO comparison. See `KEPLER_SLIP_WALL_NATIVE_2026-10-29/README.md`.

## Other completed acceleration evidence

A100 job **439115** finished the seven native cases but its aggregate exact geometry gate failed. As on Ada, control40 versus vertex40 passes all four geometry/RHS comparisons; halo20 differs in tiny face endpoint rounding while all final native physical fields and traces match. The failed halo20 gate is retained; no acceleration was promoted on that basis. The independent Ada halo40 reuse+alias gain remains a hardware/interval-specific result, not a demonstrated A100 production gain.

Artifacts: `KEPLER_Q1_T3_SUPPORT_COMPARISON_2026-10-29.json`, `KEPLER_Q1_T3_SMOOTH_ANNULAR_MASS_2026-10-29.json`, `KEPLER_Q1_T3_DENSITY_SUPPORT_2026-10-29.png`, native strict/radial gate receipts and `KEPLER_VERTEX_HALO_A100_RESULT_2026-10-29.json`. The reproducible readers are `audit_kepler_t3_support_20261029.py` and `audit_kepler_smooth_annular_mass_20261029.py`; four small geometric tests pass, including comparison against explicit periodic images and rejection of invalid image certificates.

The long-term objective remains unachieved: t>6 stability and matched GIZMO t=120/t=600 density fields, plus the explicitly interpreted approximately 1000-orbit statement, are still required.
