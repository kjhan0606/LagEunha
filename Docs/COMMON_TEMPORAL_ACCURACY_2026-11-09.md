# Current endpoint operator: temporal and shock controls — 2026-11-09

These are new controls of the **current common endpoint-volume coupling**, not repetitions of the archived pre-inverse CFL benchmark. No physical force, limiter, tolerance, weight cap, or live reference trajectory was changed. Previous goal turn made progress by completing and assessing KH128/Sedov64, then qualifying the whole KH64 C=.8 history. This turn adds same-start temporal evidence, area-integrated comparison, and full Sedov timestep controls. The full 1D/2D/3D accuracy/performance objective remains incomplete.

## Same-state KH temporal control

`common_endpoint_temporal_20261109_v1` restarts the identical accepted 8192-cell t=1.2 state from `common_closed_volume_20261109_v1/full_kh64_cfl08/native_t12/khout.000395.dat`, hash ff2d4e100cd51f65576d29db7afcb10708c7528aae65830b3aafb2e20cd54478. All four cases use binary c3554c696ec0b2044ac34ba20f54c4dfe7c920bf142eb6fd3e67a91f799a9ffb, endpoint inverse, cap .001, and identical initial RHS (all stored fields bitwise equal). Only Courant differs. Each reaches exactly1.22, retains every mass, positive U/V, volume error <=6.79e-12 and total energy to roundoff.

| Courant | Accepted steps | Retries | Median dt | Shared-host process seconds |
| --- | ---: | ---: | ---: | ---: |
| .8 | 28 | 13 | .000722580 | 8.777 |
| .4 | 56 | 30 | .000361452 | 14.394 |
| .2 | 111 | 52 | .000180811 | 26.943 |
| .1 | 222 | 83 | .0000904144 | 49.317 |

Relative to the smaller-step C=.1 result (not an exact solution), same-material-ID vy L2 differences are .02742%, .02073%, .01536% for .8/.4/.2. Global RMS differences are much smaller and are not substituted for field differences. No asymptotic time-order claim follows from this finite interval and the limited reconstructed operator.

## Polygon-intersection field norms

The pixel raster was sensitive to subpixel interface displacements at these small differences. `common_endpoint_temporal_20261109_v1/polygon_overlap/audit.py` therefore reuses the existing RT half-plane clipping and omitted-site certificate, extended to both periodic axes, and integrates the saved piecewise-constant fields over polygon intersections with Shapely/GEOS. It does not reconstruct or smooth fluid values. Periodic uniform-cell volume, self-comparison, and a translated discontinuity with known analytic squared L2 are checked. Every actual reconstructed cell is compared with its native saved volume; every fragment's intersection coverage is checked.

Relative C=.8/.4/.2 versus C=.1 spatial vy L2 differences are **.05354%/.04115%/.03277%**, and density L1 differences are2.52291e-6/1.09055e-6/5.83144e-7. These include displacement of cell boundaries and therefore differ from material-ID norms. They integrate numerical piecewise-constant fields with floating-point geometry, not exact continuum solutions. Geometry/coverage certificates and input/source hashes are saved in the report. This short interval cannot establish the cause of the entire long-history difference.

The same area-integration applied to the prior/current C=.8 full histories at1.35 gives density relativeL1 **.000802409** and vy relativeL2 **.0499461**, with total overlap area2, fragment coverage discrepancy9.17e-13, and native-volume differences below2.87e-14. Thus the previously reported ~4.96% raster vy difference is real for these saved numerical fields (~5.00% by polygon integration), not merely raster noise. Evidence: `common_closed_volume_20261109_v1/full_kh64_cfl08/polygon_overlap_t135/comparison.json`. The two builds/histories still differ, so do not label this pure temporal error.

## Full Sedov64 controls

`common_sedov_timestep_20261109_v1` retains the same frozen portable binary, finite source width .0625, all initial fields, endpoint inverse and final time .1 as the completed baseline. The effective baseline is **Courant=.4**, cap .000125, verified from params; do not confuse it with older C=.2 shock screens. First increase only the cap to .0005; then change Courant to .8 at that same cap. All initial RHS body values remain bitwise equal; only the diagnostic budget header changes. Both candidates reach the target with fixed masses, positive U/V, unchanged volume acceptance, and relative energy drift <4.53e-13.

| Setting | Steps/retries | Median dt | Density reference L1 | Pressure reference L1 |
| --- | ---: | ---: | ---: | ---: |
| C=.4, cap .000125 (baseline) | 827/20 | recorded in baseline log | .10789767 | .08227773 |
| C=.4, cap .0005 | 429/89 | .000253733 | .10759119 | .08244033 |
| C=.8, cap .0005 | 315/65 | .000330126 | .10845307 | .08300345 |

Reference numbers integrate the ideal point-source solution over actual cells at quadrature32. That continuum reference is not identical to the finite square source. At C=.8 the pressure difference increases by .00072572 absolute (~.88% of the baseline diagnostic error); density difference increases by .00055540 absolute. Hence fewer steps do not imply increased accuracy. Measured process seconds154.552/91.449/66.672 are shared-device costs, not an isolated timing benchmark or GIZMO comparison.

At density threshold2.5, baseline and cap-only cases have all256 angular rays, with radius RMS1.1104% and1.0770%. C=.8 has254/256 rays, so a complete fourfold coefficient is intentionally not reported. At threshold2.0 all256 rays exist in all cases. Missing threshold crossings indicate locally reduced density peaks, not proof of a literal gap in the shock surface. Morphology is not certified by conservation or a mean radial norm alone.

## Decision and next evidence

Keep enlarged-timestep settings as candidates; do not alter current cylinder/RT histories or globally adopt a CFL from these controls. The new endpoint inverse remains compatible with these larger steps in the tested KH/shock states, with measured accuracy tradeoffs. Current-node timing at matched accuracy remains unproven, as do remaining 2D morphology/refinement, MPI>1, and eventual3D. The next Sedov spatial comparison should use the already prepared fixed-width128 IC and the baseline timestep scaling, preserving the physical source rather than shrinking it with particle count. That separates source/lattice anisotropy and resolution effects without adding a force or smoothing patch. No3D evolution or second-paper drafting was started.


## Fixed-source Sedov128 refinement launched

`common_sedov128_20261109_v1` now evolves the prepared128² fixed-width source from0 to.1, with the same portable binary and C=.4 endpoint method as64²; dtcap halves to.0000625 with grid spacing. Actual extracted mass/U/V aggregate four fine cells to each coarse cell within roundoff; velocities and W are zero in both. Checkpoint allocation is at most4 including the copied input (~45.6MB), dump interval10, CPU12/13 and the available third Ada slot. NativePID300606 and supervisor300582 were verified live with accepted evolution beyondt.00136. The finite runner checks unchanged invariants and then invokes the existing cell-integral/angular analysis with only its resolution assertion changed to16384. This is a new spatial refinement, not a completed result, and no3D is launched.

## Current KH128 full-history temporal control

After Sedov128 evolution completed successfully (CPU shape analysis continues), launched `common_closed_volume_20261109_v1/full_kh128_cfl08` on the freed Ada evolution slot, CPU12,13. This is a new 128×256 temporal control, not a repeat of completed 64×128 experiments. Same canonical binary SHA c3554c696ec0b2044ac34ba20f54c4dfe7c920bf142eb6fd3e67a91f799a9ffb, pristine source, physical operator, endpoint tolerance1e−11, reconstruction, EOS and dtcap.001 as completed current KH128C=.2. Parameter file differs only in Courant_hydro .2→.8. Memo/LU/symbolic reuse remain off. Copied source all13 fields and full initial RHS (including header) are bitwise equal.

Targets .3,.6,.9,1.2,1.35. Maximum15 checkpoints including copied inputs,342099660bytes; measuredfree149537289142272bytes. Effective params and environment recorded per segment; dump interval10 exceeds target; original live cylinder/RT unchanged. Initial accepted dt~.000382636, four times original KH128dt. Invariants are enforced at each target; completed physical/timestep comparison still required. No global CFL adoption or matched-accuracy≤2× performance claim. Launch supervisor322176, native322181, execsession38387.

### KH128 first matched endpoint, t=.3

Current C=.8 completed785steps/105retries versus baseline C=.2 3137steps/137retries. Transverse KE relative difference−3.65683e−5 and fundamental-amplitude difference+1.66621e−5. This is an early temporal comparison, not an accuracy pass. Summed accepted dt=.3 in both histories; max reported time-increment discrepancy<2.8e−17, so the observed generator differences are not a mismatched output timestamp. Shared wall354.331s versus1191.618s is not a controlled performance benchmark. Full1.35history remains running. Evidence: `full_kh128_cfl08/early_comparison.json`.

## KH128 spatial fields and independent early reference

The temporal-analysis worker (PID344824, execsession36053) uses the previously qualified periodic power-cell clipping and omitted-site certificate. It integrates piecewise-constant fields over polygon intersections, including periodic fragments. Analytic translated-discontinuity and exact-self fixtures pass again. This is a CPU postprocessor on16,17, not another GPU evolution. It evaluates completed .3/.6/.9/1.2/1.35 endpoints as available, retaining the native calculation unchanged.

| Matched time | Density relative L1, C=.8 vs C=.2 | Transverse velocity relative L2 | Low y-momentum modes relative L2 |
|---|---:|---:|---:|
| 0.3 | 0.004807% | 5.4819% | 0.2828% |
| 0.6 | 0.006428% | 2.0705% | 0.3617% |
| 0.9 | 0.009790% | 1.1092% | 0.2745% |

At.3, the5.48% full-field difference contrasts with the0.0017% fundamental-amplitude difference previously reported. Both are correct but measure different things. Same-ID vyL2 is5.475%, and absolute velocity RMS difference is1.982e−4;94.46% of its squared contribution comes from material initially within.1 of either interface. At.6, same-ID vyL2 is2.007%, absolute1.714e−4, with87.31% in that initial band. Geometry relocation alone is therefore insufficient to explain the field difference. These band and same-ID diagnostics are not Eulerian norms and do not by themselves prove numerical noise versus physical small-scale structure. Total thermal-energy difference is approximately−1.4e−8 relative at both early times; this field discrepancy must not be conflated with gross global heating.

To assess whether C=.8 is actually less accurate, reused the existing linear Euler reference (`kh_linear_reference_20261106_v1/solve_density.py`) at previously unsaved .3/.6 times; only its time array/end time changed. It remains a CPU validation calculation, not the production solver.512/1024 transverse reference energy changes by<1.43e−8, and1024 agrees with the independent compiled-C reference within6.12e−10 relative. The same existing first-order conserved-reference formulas were integrated on each native mesh with3/5-point triangle quadrature. Reference field changes are<4.51e−9 absolute, and cell volumes match within~7e−14.

| Time | C=.2 vyL2 vs cell-integrated linear reference | C=.8 vyL2 vs same reference |
|---|---:|---:|
| 0.3 | 4.20499% | 4.30527% |
| 0.6 | 3.04828% | 3.08029% |

The large pairwise field difference is not the increase in error against the independent reference: early reference discrepancies are similar for both timesteps, though C=.8 is slightly larger. Neither reference comparison is an exact nonlinear error: finite-amplitude effects, generator-sampled initial velocity, and first-order perturbative conserved fields limit the interpretation. Do not adopt C=.8 solely from growth amplitude or reject it solely from pairwise L2; finish the full trajectory, shape comparison, and same-accuracy performance check. No limiter, viscosity, weight rule, or physical coefficient changed. C=.8 remains a candidate, not a global default.

Evidence: `full_kh128_cfl08/temporal_analysis/{run.py,comparison.json,status.json}`, `local_temporal_difference/{run.py,comparison.json}`, and `linear_cell_average/{run.py,reference.py,comparison.json}`. The four early cell-reference analyses completed successfully; the full-history analysis worker remains active.


## KH128 C=.8 final temporal assessment, t=1.35

The full 128×256 candidate and its independent CPU postprocessor completed. The final actual-power-cell density figure was inspected: both settings retain the same large-scale interface wave; this is not yet a developed KH roll-up. No force, limiter, viscosity, weight equation, tolerance, or physical coefficient changed. Only the Courant coefficient differed from the current C=.2 reference.

| Time | Density relative L1 | Transverse-velocity relative L2 | Low y-momentum modes relative L2 |
|---|---:|---:|---:|
| 0.3 | 0.004807% | 5.48191% | 0.28284% |
| 0.6 | 0.006428% | 2.07047% | 0.36175% |
| 0.9 | 0.009790% | 1.10923% | 0.27449% |
| 1.2 | 0.021704% | 2.25757% | 0.25007% |
| 1.35 | 0.042676% | 2.90113% | 0.31942% |

These are temporal differences from C=.2, not errors against an exact solution. The final density difference is small, but the 2.90% transverse-field difference must not be replaced by the smaller integrated-mode statistic. Certified final geometry covers area2, native-volume agreement is better than7.2e−14, and omitted-site margins remain positive for both meshes.

Against synchronized GIZMO, low y-momentum mode discrepancies are 15.77049% for C=.2 and 15.79942% for C=.8. Thus the timestep enlargement does not resolve the main cross-method growth discrepancy. GIZMO is not an exact solution, and its post-initialization pressure offset (~0.0581%) and velocity change remain documented. The earlier cell-integrated linear-reference checks support comparable early accuracy, not a nonlinear late-time accuracy certificate.

Every stored endpoint retains exactly the original cell masses and positive U/V. Maximum accepted endpoint-volume mismatch is 9.9997011e-12; maximum absolute relative energy drift is 6.7994476e-16. These do not by themselves certify morphology.

Decision: retain C=.8 as an efficient, full-history-tested KH candidate; do not change global defaults or the live cylinder/RT histories. The next performance certificate requires current-code native/GIZMO timing on the same execution resources together with the existing accuracy assessment. Remaining 2D shape/refinement and boundary qualifications still precede 3D. No second-paper drafting was started.

Evidence: `full_kh128_cfl08/results.json`, `temporal_analysis/comparison.json`, `temporal_analysis/t135/density.{png,pdf}`, and `cost_audit/report.json`.
