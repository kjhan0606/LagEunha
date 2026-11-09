# Common-method Sedov validation, 2026-11-09

The preceding status turn was a verified wait: the fine cylinder native process was live. This batch advances the common-method validation without changing or restarting the cylinder, KH, or RT runs. No 3D calculation is launched.

## Initial-condition audit

The existing C writer was reused. Newly generated 32 and 64 square-grid ICs, with both the original four-cell source and fixed physical source width, reproduce all four archived checkpoints byte for byte. The 128 ICs were generated but have not evolved. At zero weights, the verified uniform square lattice has exact rectangular Voronoi cells; mass, volume, and thermal energy are cell integrals.

The default source deposits unit excess energy in the central four cells: width 2/N, hence .0625 at 32 and .03125 at 64. It changes the physical initial condition with resolution. The already existing fixed-source option retains width .0625: 4, 16, and 64 hot cells at 32, 64, and 128. Restricting the fine cells onto the coarse grid reproduces mass, thermal energy, and volume exactly at both refinement levels.

Exact rectangular integrals of excess-energy-weighted Re[(x+iy)^4] give -w^4/60; the r^4 moment is 7w^4/180. Their ratio is -3/7 for the finite square source. This is a fourth moment, **not** an angular Fourier coefficient. It shows initial fourfold structure, not proof that all final anisotropy originates there. Particle-point moments would describe a different field and must not replace cell integrals. The initial source centroid is (.5,.5), density is uniform, initial velocity and weights zero, excess energy one, background pressure 1e-6, and gamma float32(5/3).

## Matched full-history comparison

The portable qualified Ada binary is frozen with SHA256 11a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e. Both 64-squared cases start at t=0 from the same archived fixed-width IC (SHA256 a7735763377c6046af3986fa28ed6f75ff7e823ea7b287cbb3dac51dbab79925), evolve to .1, and retain the same reconstruction, boundary, EOS, and timestep cap .000125. Only the endpoint volume-coordinate option differs. The control reports its volume mismatch; only the endpoint case must meet the existing 1e-11 volume criterion. Both require exact retained masses, positive finite states, periodic volume closure, and energy conservation. No source-specific hydro term, limiter, or entropy repair is added.

A finite supervisor waits for the already running full KH128 calculation to finish successfully and verifies the recorded native process has exited before using its Ada slot and CPU12,13. It does not stop or restart existing calculations. Actual status, PID, environment, command, lineage, hashes, and results live under the batch directory. A wait is not an executed Sedov time integration. The two runs execute sequentially; no fourth simultaneous GPU calculation is authorized by this runner.

The intended output is two final snapshots. Conservative allowance is six new snapshots plus two copied inputs, about 22.8 MB. Dump period 10 exceeds target .1; the runner rejects more than four checkpoints per case. Available storage at preparation exceeded 149 TB. Effective namelist is common_sedov_20261109_v1/params.dat and copied to each run directory.

## Remaining assessment

Measure the final shock contour and angular variation, pressure/density against actual-cell-integrated radial profiles, and control-versus-endpoint differences. The point-explosion Sedov solution is a diagnostic, not the exact finite-square-source solution. Source refinement and mesh refinement are separate questions. Neither the IC audit nor invariant completion proves morphology accuracy, dimensional generality, or the required matched-accuracy GIZMO speed. The first manuscript is not updated with unqualified results, and the second manuscript remains deferred.

Artifacts: Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_sedov_20261109_v1/{initial_audit.json,plan.json,prepare.py,run.py,status.json}.

## Completed postprocessing qualification

The saved-state assessment now reuses the existing C periodic power-cell reconstruction and analytic Sedov profile. On the archived fixed-source 64-squared predecessor, independent refactoring reproduces stored cell-integrated density/pressure reference metrics at quadrature orders 16 and 32 within 1e-12. This is postprocessing qualification, not a new hydro result. The order-32 density and pressure diagnostics are .1089178866 and .0844744220 for that OLD state.

Ray/polygon density-threshold measurements were added with three reported thresholds, not a selected best-looking one. At threshold 2.5 all 256 rays intersect the old shock shell; radial RMS variation is 1.191% and fourfold contour amplitude 1.281%. At threshold 3, only 222 rays intersect, explicitly recorded rather than interpolated. These contour diagnostics are not exact shock positions or proof of the origin of anisotropy.

The finite CPU assessment worker waits for the two matched new trajectories to finish and pass their invariants, then checks source/binary/environment equivalence before comparing actual polygons and cell-average radial references. Its live handle and results are under analysis_pipeline/. Evolution was still waiting for the existing KH GPU slot at launch. No extra GPU calculation was started by this analysis worker.

## Paired full-history completion and shape assessment — 2026-11-09

Both 64² finite-source cases completed exactly t=.1 from identical initial fields and the same frozen portable binary. Control: 817 steps/17 retries; endpoint inverse: 827/20. Both preserve every mass, positive U/V, total volume1 and energy relative drift below7.1e-13. The largest time-update target/geometric volume mismatch falls from6.33114e-4 to9.67118e-12. This is distinct from the independent final polygon/native-volume cross-check (~1e-13 for both).

The CPU analysis completed, with figures visually inspected at `assessment/density.png` and `radial.png`. Actual-cell integrated point-source reference relative L1 at quadrature32: density .10912089→.10789767 and pressure .08279679→.08227773. At density contour2.5 all256 rays are present; radius angular RMS1.2074%→1.1104%, fourfold amplitude1.3108%→1.1922%, mean radius .3680485→.3678224 (ideal point-source radius .3648595). The improvement in shape is modest, not elimination of anisotropy. Density contour3 misses26/28 rays, explicitly excluded from complete angular-harmonic claims.

The finite square energy-deposition region differs physically from an ideal point explosion, and the initial lattice/source anisotropy remains visible in both final maps. Thus the integrated reference discrepancies are useful diagnostics, not exact errors for identical continuum data. No force/viscosity/weight cap was added to hide the pattern. The paired test shows geometry consistency can be greatly improved while leaving much of the shape error intact. Resolution/source-geometry comparisons remain required; 2D is not complete and no3D evolution is justified by this result.

## Fixed-source 128² refinement completed

`common_sedov128_20261109_v1` reached t=.1 in1682steps/29retries,647.603s on shared Ada. Exact cell masses retained; positive U/V, volume1, energy relative−5.378735e−13, max evolution volume residual9.751304e−12. Same frozen portable binary and sourcewidth.0625 as64², C=.4 and dtcap scaled .000125→.0000625. Fine IC aggregates exactly into coarse m/U/V. No new physical operator or smoothing.

Actual convex polygon integration independently matches stored volume within1.341223e−13. With32-point triangle quadrature, point-source diagnostic L1 density10.78977%→8.73850%, pressure8.22777%→8.15257%. Fine16→32 quadrature changes density/pressure diagnostics by9.52e−6/9.74e−6 absolute, smaller than resolution change. Ideal point source differs from the fixed finite square source, so these are not exact finite-source solution errors.

At rho=2 threshold (all256rays present at both resolutions), relative angular RMS1.19662%→1.17417%, but fourfold amplitude1.11710%→1.52187%. At rho=2.5, fine has254/256validrays, so do not report a complete fourfold coefficient or mistake the two missing threshold crossings for literal holes. Thus sharper mean shock/density profiles do not establish shape convergence: residual directional structure persists and the fourfold component increases. Both density and radial figures were visually inspected; shell overshoot and interior directional troughs remain. No shape/physics pass.

Evidence: `common_sedov128_20261109_v1/comparison.json`, `assessment/report.json`, `coarsening_check.json`, and `figures/{density,radial}.{png,pdf}`. Plot uses actual polygons with shared rho range0–5; dashed circle is ideal point-source radius. Retain both data sets. Next source-shape/mesh diagnosis should reuse prior IC controls where available, not add a shock-specific damping patch. Full KH128 temporal control now occupies the freed GPU slot.
