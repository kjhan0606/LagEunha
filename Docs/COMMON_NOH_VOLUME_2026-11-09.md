# Common geometric-volume coupling for the existing Noh exterior

The preceding goal turn made progress through actual high-resolution saved-state reuse comparisons and the planar Sod evidence. This turn removes the Noh exception from the already-qualified common volume audit/inverse. It does not introduce a new shock model, boundary state, entropy repair or tolerance. The existing upstream analytical exterior and owned-cell membership are retained. Only single-rank native2D is qualified; no3D or general MPI claim follows.

## Operator and boundary contract

The exterior RHS bridge optionally returns the realized volume rate that already enters `Udot=-P*Vdot+H`. The native Noh adapter samples that exact rate into the existing SSP stage audit. With the audit disabled the additional output is null; with audit enabled but inverse disabled all13 saved state fields reproduce the old executable bitwise in both short fixtures.

For fixed generator positions during the Newton weight correction, each resolved power face contributes `A/(2d)*(delta_W_i-delta_W_j)` to the owned-cell volume change. Owned/owned faces use the neighbor correction. An exterior sample is prescribed, with fixed W=0 and delta_W=0, so it supplies an anchor. The inverse reuses the SAME exterior resolved-face builder as the pressure-work calculation, including existing paired tiny-face handling and geometric certificates. Exterior samples are never added to owned mass or promoted to evolving fluid variables. No periodic aliases or closed-domain gauge projection are applied to this open-boundary problem.

The already qualified anchored solver retains its1e-10 relative raw-equation residual certificate. Rebuilt actual GPU geometry must meet the unchanged1e-11 relative target-volume criterion before acceptance, with the existing eight-iteration cap. The same endpoint or predictor+endpoint choice is available as for the other supported2D flows. No force, momentum, mass, internal energy, prescribed sample coordinate or external weight is modified by the volume solve.

## Qualification

Ten short native runs use the same saved32²t=.1 and64²t=.4 states. Reference/off/audit/endpoint/stages all reach the same requested endpoints, without retries, with fixed masses and positive U/V. Off and audit match reference bitwise. The64² interval contains a negligible final exact-stop remainder (9 accepted steps, versus8 nominal steps); this is retained in the record.

Six extended comparisons use the same starts and dt=.000125 over a.02 interval. All finish with zero retries and no failed linear corrections:

| Grid | Mode | End time | Steps | Retries | Maximum volume error | Linear corrections | Maximum step boundary-energy residual |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| n32 | audit | 0.120 | 160 | 0 | 6.3852e-09 | 0 | 3.0575e-16 |
| n32 | endpoint | 0.120 | 160 | 0 | 9.0522e-16 | 160 | 3.0466e-16 |
| n32 | stages | 0.120 | 160 | 0 | 1.4776e-12 | 320 | 3.2179e-16 |
| n64 | audit | 0.420 | 161 | 0 | 1.0154e-05 | 0 | 6.0433e-16 |
| n64 | endpoint | 0.420 | 161 | 0 | 8.8728e-12 | 160 | 6.0802e-16 |
| n64 | stages | 0.420 | 161 | 0 | 1.0310e-13 | 640 | 6.1149e-16 |

The64² endpoint includes a2.22e-15 final remainder step. The much smaller volume mismatch establishes coupling consistency, not better shock radius, density plateau or central heating. The `stages` option costs more corrections; these short tests do not establish a need for it in the full trajectory. Boundary energy residual is an absolute per-step balance, not an accumulated error or complete entropy diagnostic.

KH, RT, Sedov and fixed-solid cylinder endpoint regressions reproduce all13 saved fields bitwise. The first regression launch failed before any hydro step due to CUDA out-of-memory while the extended Noh test also used the shared GPU; it was preserved in `regression/terminal.json`. Once the extended calculation terminated, `regression_sequential` passed all four cases. The original cylinder/KH processes were not stopped or restarted. Future additional GPU qualifications must run sequentially alongside those two persistent jobs.

The four-file change is now adopted in canonical `Exam/NativeU64`: exterior RHS output, Noh sampling, volume inverse boundary mapping and SSP activation. `adoption.json` records before/after hashes and132 local compiled dependencies verified against the adopted source. The qualified executable47ef8b9704c4970a24760461f07d37bd2f4658518d18660c760a6b46f39232ea remains immutable; existing production binaries are unchanged.

## Whole-history check launched, not yet a physics pass

`full_noh64` starts from the original4096-cell IC with SHA256a8419079e2c8dfb37957acf232d747a8b91ef561f89d70f33d1b2390a2b56053. It advances to.2,.4,.6 through verified successful checkpoints. Endpoint inversion is on; symbolic/previous-LU reuse is off to isolate the method. The old IC, reconstruction, pressure/energy flux, exterior state and CFL settings are retained. Three final-output checkpoints plus copied restart inputs total about17.1MB. CPU14/15 and the shared Ada are used. A finite supervisor records actual process handles and stops on a failed invariant; launching is not completion. Compare the final density profile and angular structure to previous native evidence and the strong-shock Noh reference before claiming accuracy. Current finite initial pressure and cell averaging must be acknowledged. Refinement, general MPI,3D and matched-accuracy GIZMO runtime remain outstanding. The second-paper draft still waits for the first astronomy draft.

Evidence root: `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_noh_volume_20261109_v1`. Exact patch/build, saved-state lineage, effective environments, logs, short/extended outputs, regression failure and recovery, adoption manifest and full-history plan are retained.


## Matched whole-history control and early spatial evidence

The first new-method checkpointt=.2 has completed and its4096 actual power cells were independently reconstructed using the unchanged prescribed exterior trajectories. Omitted-site distance/weight certificates are positive; all cells have positive area and convexity, and the maximum area/native-volume difference is1.558e-13. This is saved-state geometric verification, not an export of the native face array. Density images fill these polygons with their saved cell averages, without interpolation/smoothing. Polygon edge antialiasing is disabled to avoid white raster seams that could look like physical gaps. The red dashed circle marks the cold-reference shock radius.

For gamma=float32(5/3), the cold2D Noh reference has shock speed(gamma-1)/2, density behind the shock((gamma+1)/(gamma-1))^2 and pressure(gamma-1)*rho_post/2. This follows from radial upstream compression and the strong-shock jump conditions for unit inward speed and initial density. The simulation uses finite initial pressure1e-6 and the existing finite-pressure exterior approximation, so this is explicitly a cold-limit reference rather than an exact solution of the finite-pressure problem.

At t=.2 the cold shock radius is.066666663. All176 cells wholly inside that circle have volume-weighted mean density14.91195 (reference16.0000014), range9.51076--17.59610, and mean pressure5.07935 (reference5.3333335). The102 cells wholly inside.75 times the shock radius have mean density14.91837. These departures remain despite geometric-volume agreement. They are not explained away by energy conservation, nor attributed to a new physical entropy defect from this comparison alone.

A separately defined shape diagnostic intersects256 radial rays with the actual high-density polygons at the midpoint between the cold pre/postshock densities (rho_threshold10.0000008). All rays find a contour, none reaches the2*shock-radius search limit. Mean radius.06588664 is1.17% below the cold shock radius; angular RMS/mean is2.42%, and the fourfold harmonic amplitude is0.0496%. These measure the finite-cell threshold contour, not an exact shock location or a cell-average density error norm. Radial plot dots are values at volume centroids; cells crossing the shock are not incorrectly scored against point exact values.

The old Noh64 record used a transferred intermediate state and an earlier executable, so it cannot isolate the effect of endpoint volume correction over an entire history. A new matched control is queued AFTER the currently running full_noh64 completes. It uses the identical frozen binary, original IC hash, output times .2/.4/.6 and timestep settings, but leaves endpoint inversion off while retaining the volume audit. The control intentionally does not require a small geometric-volume defect: measuring that difference is its purpose; positive U/V, fixed masses and successful endpoints remain required. It adds three checkpoints (about17.1MB including copied restarts) and no concurrent extra GPU allocation. `matched_control_pipeline` verifies the existing supervisor/terminal results and stops if the dependency fails. It never restarts the existing calculation.

A finite CPU-only `matched_analysis_pipeline` waits for successful control completion, then reconstructs paired polygons and generates density/profile comparisons at all three times. Same IC/binary hashes and endpoint times are checked before comparison. Planned final artifacts are `comparison_t06/density.png`, `radial.png`, and `report.json`; they do not exist as completed comparisons yet. Generation is not scientific acceptance. Sources, current handles and stages are in the two pipeline directories; `analyze_history.py` is the frozen analysis script. Current preview/report lives under `preview_t02` and cumulative energy balance explicitly subtracts all exterior work since IC. No hydro operator or tolerance changed during this analysis.


## Completed matched Noh histories: consistency gain without morphology gain

Both 64² histories reached t=0.6 from the same pristine IC and frozen executable. All 4,802 accepted steps have identical accepted times and zero retries. Masses are unchanged; U and V remain positive. Endpoint inversion is the only changed method option. This supersedes the earlier pending status above.

The maximum stepwise relative mismatch between the pressure-work target volume and rebuilt geometric volume falls from 4.7971e-3 to 9.7721e-12. The largest control defect occurs at t=0.02475, near the earliest central compression, rather than at late time. This diagnostic differs from independently reconstructing the *saved* power-cell volumes: the latter agree with saved volumes at about 2e-13 in both runs and do not test pressure-work consistency.

However, final morphology is almost unchanged. The cold-reference radius is 0.199999988. The density-threshold contour has mean radius 0.202401803 without inversion and 0.202402948 with inversion (both about +1.20%). Angular RMS/mean is about 1.472% in both. Fully postshock volume-weighted mean density is 15.499817 versus 15.499575, compared with the cold reference 16.000001; mean pressure is 5.191204 versus 5.191111, compared with 5.333333. Directional internal bands and a central density deficit remain visible. These are not removed by enforcing endpoint volume consistency. The new-minus-control, mass-weighted density difference at matched material IDs is only 0.014744% at t=0.6; this is not an exact-solution spatial error norm.

Cumulative energy balance after subtracting exterior work is 2.9320e-13 (control) and 2.9404e-14 (endpoint), both very small. This alone cannot certify morphology or entropy. The original finite-pressure/exterior and cold-reference caveats still apply.

The evidence therefore supports retaining geometric/pressure-work consistency as an explicit invariant, but does NOT support claiming it resolves Noh's main thermodynamic or shape errors. Do not add a Noh-specific correction on this evidence. Next physics diagnosis must distinguish initial lattice/central radial-velocity representation, reconstruction/dissipation, and boundary effects, using the existing common operator. Full RT/Sedov morphology and refinement remain incomplete; 3D is still deferred until the 2D suite is qualified.

Evidence: `matched_assessment/report.json` and `comparison_t02`, `comparison_t04`, `comparison_t06` (density, radial plots, reports). Both final images were inspected. Plots are diagnostic, not declared publication-ready. No force, EOS, reconstruction, dissipation, timestep tolerance or IC was modified for this assessment.

Observed shared-GPU wall totals were 552.80 s without inversion and 705.61 s with inversion, ratio 1.276. These sequential, shared-hardware timings are not an isolated performance comparison or a GIZMO benchmark.
