# Sedov source-orientation control — 2026-11-09

Previous turn completed the same-width64²/128² current-method histories and found that density-reference mismatch decreased while fourfold directional structure did not. These are new causal controls of that remaining2D shape issue, not a new hydro model, smoothing operation, or second-manuscript campaign. First-paper validation remains the priority;3D evolution remains deferred.

## Minimal IC change

`Exam/NativeU64/unified_initial2d.c` now accepts optional `GFS_SEDOV_SOURCE_ANGLE_DEG` for the Sedov source only. A convex half-plane intersection integrates the rotated physical square into each original rectangular cell. The mesh, generator locations, masses, volumes, velocities, weights, total injected energy, source width, EOS and solver are retained. Default/zero angle uses exactly the original arithmetic. Multiples of90 degrees are reduced to zero for the square. Invalid nonfinite/malformed angles fail before writing a checkpoint. The exact before/after sources and hashes are preserved in the audit directory. Running binaries were not rebuilt or changed.

The source remains a finite square of width1/16 and total excess thermal energy1. Angle45 is a diagnostic input, not a viscosity/weight/limiter fitting coefficient. Rotating the physical source relative to a fixed lattice may also change its finite-cell projection; that difference is measured, not hidden.

## Completed IC verification

- Recompiled before/after writers produce bitwise identical default KH, Sedov, Gresho, periodic shear, periodic Sod and Noh checkpoints. These are six IC parity tests, not six new hydrodynamic evolutions.
- Angle0 at32/64/128 and angle90 at64 reproduce the archived fixed-source checkpoint byte for byte.
- Independent Shapely polygon intersections at32/64/128, including22.5 and45degrees, reproduce per-cell deposited energy to at most2.26e−16 absolute.
- Coarsening rotated45-degree fine cells into coarse cells reproduces mass/volume exactly and thermal energy within2.78e−17 absolute. Total energy and centered centroid are retained to roundoff.
- Projected fourth-moment ratio H4/R4 at45degrees is+.29492 (64²),+.39165 (128²), approaching continuous+.428571; aligned source is−.428571 exactly. This ratio is not an angular Fourier coefficient. The projected second radial moment is13.42% high at64 and3.22% high at128 relative to the continuous rotated square. Hence the64 result alone cannot cleanly separate source rotation from projection width.

## Finite evolution plan and actual state

Runner `common_sedov_orientation_20261109_v1/run.py`, supervisor335996 / execsession79502, is waiting for the existing KH128C=.8 full history to release the third Ada hydro slot. Waiting is not a Sedov time evolution. It checks the actual dependency executable/handle and successful terminal results; it neither restarts nor stops KH, RT, or cylinder. A six-hour wait bound fails explicitly if unmet. Syntax jobs456829/456838 remain pending QOSGrpGRES/Dependency, so no duplicate H200 submission was made.

Once free, run45-degree64² then128² to t=.1 using the same portable frozen hydro binary as aligned source controls (SHA11a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e), C=.4, caps.000125/.0000625. CPU12,13, GPU0; CPU analysis20,21. At most4checkpoints including input per case, estimated57017248bytes total; measuredfree149534739005440bytes. Dump interval10 exceeds target. Positive state, fixed mass, volume tolerance1e−11 and global energy checks are unchanged.

Existing polygon-based Sedov analysis follows each case. Compare full density maps plus **complex** fourth angular contour coefficient, not merely its magnitude: a rotated-source contribution should change phase under source rotation, while a grid-aligned contribution need not. Mixed/nonlinear contributions and missing density-threshold rays must remain explicit. No finite-source exact-solution, morphology-pass, or matched-accuracy GIZMO speed claim follows from this setup alone.

Artifacts: `common_sedov_orientation_20261109_v1/{initial_audit.json,source_change.json,prepare.py,plan.json,run.py,status.json}`. No commit or push performed.


## Completed64 source-rotation comparison

The KH128 job completed, releasing the third GPU slot. The finite Sedov runner then completed64² through t=.1 (818 accepted steps,17 retries,154.331 shared-device seconds) and started128². Native128 PID382270 was verified live; it is not yet a completed result. The original cylinder and RT calculations remain untouched.

The64 case retains exact fixed masses, positive U/V, total area1 and maximum endpoint-volume mismatch9.7996e−12. Relative global-energy drift is−1.1565e−12. The actual power-cell polygons agree with native volumes to5.39e−14. These invariant checks do not certify shock accuracy.

At the predeclared rho=2 contour, all256 rays exist in both cases. The complex fourth coefficient changes from (+.0111709843,−1.50e−9) for source0 degrees to (−.00633166918,+1.34e−8) for source45 degrees. Its principal orientation therefore changes from essentially0 degrees to−44.99997 degrees (equivalent45 modulo90). This is evidence that the source/projection orientation contributes to the fourfold deformation; it contradicts an explanation based solely on a fixed-grid orientation independent of the source. It does NOT measure a unique additive split into source and grid errors or exclude their nonlinear coupling.

The angular RMS changes1.1966%→1.0141%, but point-explosion reference mismatch becomes larger: density L1 .107898→.120830 and pressure L1 .082278→.114706. Mean threshold radius changes .370238→.370843 (ideal point-source radius .364860). The finite square and its finite-cell projection differ from the point source, so these are reference discrepancies rather than exact finite-source errors. A smaller fourth-mode amplitude alone does not justify choosing the rotated IC as a better solver result. rho=3 lacks8/256 rays in the rotated case; no complete fourth coefficient is claimed there.

Actual-cell density and angular figures were inspected. `figures64_v2` reserves empty plot margins for the inside-panel title/legend and preserves the raw contour variations; no smoothing was applied.128² will check the same source-angle contrast with a substantially smaller initial projection error. No force, viscosity, limiter, weight law, or boundary patch was introduced. Default Sedov IC remains unchanged.

Evidence: `run64/result.json`, `assessment64/report.json`, `comparison.json`, and `figures64_v2/{density.png,angular.png,comparison.json}`.


## Interrupted128 execution and recovery

The original supervisor335996 was confirmed absent while native382270 was still live. Polling its existing exec session79502 returned terminal143. Native382270 was subsequently confirmed absent too; its log stops after accepted t=.096693923890631805 and part of the following step, with no target checkpoint. Termination initiator is unknown; the log tail contains no hydro rejection/failure diagnosis. This is neither successful endpoint completion nor evidence of an algorithmic instability.

Only the t=0 checkpoint survives in the original run128 directory, so resuming near.0967 is impossible. `run128_retry1/recovery.json` records the actual available checkpoints, old-log hash, absent handles and supervisor exit; original files are preserved. Detached supervisor402999/native403019 now executes the same128 IC, frozen binary, parameters, CFL and tolerances in a separate directory. Output interval alone changes10→.05 to retain a midpoint checkpoint; at most4 checkpoints (~45.6MB) are allocated. Original and replacement initial RHS files are bitwise identical. No64 run is repeated; original RT/cylinder remain untouched.

The replacement is detached from the tool session. Its bounded finalization invokes the existing128 analysis and publishes canonical result/state files only after the unchanged target/invariant gates pass. Those files did not exist in the interrupted directory; original logs/checkpoints are never replaced. Gresho remains waiting for successful Sedov completion, not merely GPU availability. No completed128 orientation result is claimed at this point.


## Completed128 comparison and interpretation

The recovered128 evolution and its independent saved-cell analysis have completed successfully. Endpoint t=.1 was reached in1785 accepted steps with117 retries and682.845 solver seconds. Exact retained masses, finite positive states and total area1 were preserved; max endpoint-volume mismatch9.9128e−12, global energy relative change+1.8450e−13. Native endpoint checkpoint is `run128_retry1/khout.001786.dat`. Original interrupted logs/checkpoints remain intact. This is the recovery result, not a successful completion of the interrupted process. Canonical result/state copies identify the recovered checkpoint.

At128² all256 predeclared rho=2 rays exist. Rotating the square source changes the complex fourth coefficient from (+.0152181085,+.0001361144) to (−.00605976899,+7.15e−15): orientations−.1281° and−45°. Angular RMS decreases1.17417%→.60643%, fourth-mode amplitude1.52187%→.60598%. The source orientation affects the deformation at both resolutions; a source-independent fixed-grid explanation alone is insufficient. The result does not isolate a unique additive source/grid contribution, since finite-cell projection and nonlinear coupling also change.

For128² the ideal point-source density L1 discrepancy decreases8.73850%→7.81360%, pressure8.15257%→4.72908%. This trend differs from64², where rotating the source increased both reference discrepancies. Neither source angle is adopted as a tuned numerical cure. Finite-square results are not the exact point-explosion problem. The rho=2.5 aligned contour has only254/256 rays; no complete-mode comparison is claimed at that threshold.

`figures128_v2/density.png` and `angular.png` were generated from hash-verified states and actual saved convex polygons and visually inspected. Shared density range0–5, no smoothing; dashed circle is the ideal point-source radius. The density maps still show localized deficits and fine-scale structure behind the front. A smaller contour coefficient does not establish satisfactory morphology. Labels and legends lie inside panels, without covering plotted data.

`interpretation128.json` records evidence hashes, both resolutions, scope and the next causal decision: retain the common hydro and use a rotationally symmetric source control with independently verified energy deposition and explicitly matched source extent before contemplating any operator correction. No new viscosity, limiter, force or weight patch has been added. This prospective control is not yet implemented or run. Gresho has taken the released third Ada evolution slot; cylinder and RT continue unchanged.
