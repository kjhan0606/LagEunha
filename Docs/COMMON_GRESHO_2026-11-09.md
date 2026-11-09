# Current common-method Gresho qualification — 2026-11-09

Prior64²/128² full trajectories through t=1.2 exist, but predate the common endpoint-volume inverse. They remain historical evidence and are not relabelled as validation of the latest coupling. No old time interval is restarted or overwritten.

## New bounded calculation

`Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_gresho_20261109_v1` uses the pristine existing64² IC, SHA8d5e105d80735b14f38c7402a9655217fb053836bbb80e9560b8137f8fc1c9ea. Extracted initial13 fields match the previous run bitwise. Initial density is exactly1, W=0, and velocity matches the generator-sampled analytical field to roundoff. The unchanged common portable binary SHA11a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e is the same used for the current Sedov/RT tests. No vortex-specific physics or tuning is added.

Targets .002 (initial-step qualification), .3, .6, .9,1.2. Reuse C=.4 and dtcap.00025 from the historical64 test, with current endpoint tolerance1e−11. No LU/memo acceleration. Original periodic boundary, gamma float32(5/3), mass, pressure quadrature, and generator-sampled velocity remain unchanged. Requested dynamics are common SSPRK2 despite old explanatory KDK comments in the parameter template.

Supervisor392066/session54243 currently waits on the live source-rotation Sedov job; this is not yet Gresho evolution. It verifies both completed Sedov results and the absence of its actual native handle before using the same third Ada GPU slot on CPU12,13. No fourth concurrent GPU evolution. At most15 checkpoints including copied segment inputs (~42.8MB), dump interval10 exceeds all stops; free space at preparation149322641440768bytes. Effective params, environments, hashes, and lineage are written per segment. Failed invariants stop the batch without restarting a live job.

## Assessment contract

Check exact fixed masses, finite positive U/V, area1, accepted endpoint-volume mismatch≤1e−11, and energy drift≤1e−10. Separately report analytical generator-sampled velocity L1, density L1, radial RMS and outer-flow speed. Bulk generator angular momentum excludes cell-internal spin and may include periodic-boundary torque; it is not an exact angular-momentum certification.

`assess.py` additionally integrates the steady analytical velocity/pressure over actual saved convex power polygons using16/32 quadrature. The original velocity IC is point-sampled, so its initial cell-mean discrepancy must be reported. Historical/current comparison is not a clean endpoint-off/on test because historical source differed. Report quadrature variation and preserve this limitation. No3D evolution, global default change, physical accuracy pass, or GIZMO performance certificate follows merely from completion.


## Initial and historical cell-reference assessment completed

The new CPU diagnostic reconstructs all periodic power-cell polygons and checks positive convex geometry/area agreement, then integrates the existing analytical Gresho profile using triangle quadrature. At t=0, the cell-momentum relative L1 discrepancy is0.154369% even though generator-sampled velocity agrees to roundoff. This follows from the original point-sampled IC and must not be assigned to time evolution. Initial integrated pressure mismatch is4.04e−10 relative. Against the archived independent tensor128 quadrature on identical initial fields, triangle32 mean velocities agree within6.90e−7 absolute and pressure within3.57e−10. No IC has been changed.

For the historical64² t=1.2 snapshot, cell-momentum relative L1 is4.071470%, integrated pressure L1 is0.100663%, and density-volume L1 is0.103244%. Triangle16→32 changes the momentum relative L1 by8.55e−8 absolute. Geometry matches native volumes within3.93e−14. This is a new metric of the old result, not an improvement of that simulation; it must not be compared numerically as the same norm as the prior4.988% generator-sampled error.

Hydro supervisor392066 waits for Sedov. The separate finite CPU assessment worker397230/session98762 runs on CPU16,17, with status in `follow/status.json`. It waits for actual current.3 and1.2 results, certifies terminal endpoints, computes the same cell-integrated metrics, then renders actual polygon density/azimuthal-velocity comparisons. It does not launch or restart hydro. Live process absence triggers explicit failure after a grace period; six-hour wait bounds are retained. New current Gresho physics results remain pending.


## Current common-method evolution completed

All five segments reached their targets through t=1.2 using the original IC and frozen common binary. Exact cell masses, positive finite U/V, total area1 and endpoint-volume tolerance were retained; total accepted steps 4803, retries 0, solver wall time 625.861s on shared Ada. Maximum accepted volume mismatch 9.86565e-12. Final global energy change is−6.56e−16 relative. These are consistency gates, not an accuracy pass.

| t | Current generator velocity L1 | Historical generator velocity L1 | Current density L1 |
|---|---|---|---|
| 0.3 | 0.0196465459 | 0.0196465761 | 0.000271765696 |
| 0.6 | 0.0354450439 | 0.0354441904 | 0.000520019802 |
| 0.9 | 0.0430983244 | 0.0430564671 | 0.000736115981 |
| 1.2 | 0.0500050077 | 0.0498801011 | 0.00103476564 |

At t=1.2 velocity discrepancy is5.00050% versus4.98801% historically, and density discrepancy.103477% versus.103244%. The common endpoint-volume coupling preserves nearly the same vortex evolution but does not show an accuracy improvement here. The source versions differ, so this is not an isolated inverse-off/on ablation. Generator-sampled and cell-integrated norms must remain separate.

The independent current t=.3 cell-reference analysis completed: momentum-integral L1 1.670337%, pressure-integral L1 .0383454%, density L1 .0271766%; the momentum norm changes3.48e−8 absolute between quadrature16/32. Original IC cell-momentum discrepancy is.154369%. Current t=1.2 cell-reference analysis and figures are now running under the existing finite follower397230. The completed evolution has released the third GPU slot to the disk-source Sedov control; original cylinder and RT remain unchanged.


## Endpoint cell-reference and thermal budget assessment

The saved-cell t=1.2 assessment and actual-polygon figures completed and were visually inspected. Current cell-momentum L1 is4.082231%, pressure L1 .101347%, density L1 .103477%; historical values are4.071470%,.100663%,.103244%. Changing quadrature16→32 changes current momentum L1 by3.69e−8 absolute. Power-cell areas match native volumes to3.13e−14. `follow/rho.png` and `vphi.png` show closely similar earlier/current structures: a coherent rotating core remains, with density fluctuations near the transition annulus. Plot colors span actual combined extrema rather than clipping deviations. These differences do not establish improved physical accuracy from the new volume inverse alone.

`budget_audit/run.py` matches one logged thermal budget to every accepted step of all five segments, confirms time order and endpoint changes against independently extracted snapshots, and records frozen-source hashes. The frozen integrator/ledger/projected-heat definitions equal the current corresponding source files. Initial kinetic energy .0837993562 falls by .00363155321 (4.333629%) by t=1.2, while internal energy increases by the same amount to roundoff. The recorded cumulative contributions are:

| Contribution | Energy change |
|---|---:|
| Projected acoustic heat | +.003813917241 |
| Reversible pressure work | −.0001823599005 |
| Explicit compatible SSPRK2 thermal correction | −4.129586507e−9 |
| Net internal energy change | +.003631553211 |

This identifies the discrete heat term as the main positive thermal-budget contribution; the explicit time correction is approximately one million times smaller. It does not rule out time-discretization effects through the evolved trajectory, or prove that this heat term alone initiates the velocity error. Merely changing the integrator or adding an empirical damping correction is not justified by this budget.

The mass-weighted mean material-ID entropy proxy delta ln(P/rho^gamma) is+.00044999355, range−.00071200119 to+.00521286663. About72.9% of the net mean change comes from cells initially at.22≤r<.38, and14.0% from.18≤r<.22. These are cell-mean thermodynamic proxies, not exact integrals of subcell entropy; small negative cell values must not be silently clipped or automatically called a proof of a continuum entropy violation. Initial point-sampled velocity/cell-mean mismatch remains explicit.

## Current128 refinement queued behind disk-source Sedov

A detached finite supervisor455312 is waiting for the existing disk-source64²/128² calculation to finish and its actual native handle to exit. It will reuse that third Ada slot. `common_gresho128_20261109_v1/run.py` uses the pristine historical128 initial checkpoint SHA b4b5368c3d75596d2083eba4e8d24b8e9fb0c42a4fd9943645e4312d26d06df2, verifies all13 initial fields against the saved historical128 extraction, and uses the same current portable binary and endpoint-volume tolerance as current64. This fills the current-method resolution comparison; the old128 trajectory remains separate and is not relabelled.

Targets .002,.3,.6,.9,1.2 reuse the latest accepted segment checkpoint. C=.4 and cap.000125 retain the historical128 settings; no vortex-specific model change. Output interval.1 gives intermediate recovery states, at most5 checkpoints per segment including copied input,25 total (~285.1MB); free space at preparation149.286TB. Effective parameters are recorded per `tTAG/params.dat`. New128 evolution has not yet started at this report; invariant and physical/shape/thermal comparisons are still required. No3D evolution and no GIZMO speed certificate.


##128 saved-cell assessment follower

Detached CPU follower462985 uses the identical64 analysis with only cell count/root/historical128 paths changed. It processes initial and historical128 reference states before waiting for current.3 and1.2 accepted endpoints, then produces actual-polygon comparison figures. CPU16,17, about4.30GB temporary geometry workspace; measured available host RAM813.9GB before launch. It never starts or restarts hydro. Initial128 reference analysis has completed and historical reference analysis is in progress; actual current128 hydro still waits for disk-source Sedov analysis/qualification. Six-hour dependency bound and live-handle checks remain.


## Current128 t=0.3 refinement and thermal comparison

The current128 run has reached t=0.3 and is continuing to1.2. The same-time comparison is saved with source-report hashes in `gresho_glass_evidence_20261109_v1/gresho_t03.json`. No hydro terms or parameters were tuned. Both resolutions use the same common binary; the timestep cap is refined with spatial resolution, so this is not an isolated spatial-order measurement.

| Metric at t=0.3 | 64² | 128² | Coarse/fine |
|---|---:|---:|---:|
| Kinetic-energy loss relative to initial kinetic energy | 0.867691% | 0.222174% | 3.90546 |
| Integrated acoustic heat | 0.000761297881 | 0.000196850926 | 3.86738 |
| Generator velocity relative L1 | 1.96465% | 0.812150% | 2.41908 |
| Cell-integrated momentum relative L1 | 1.67034% | 0.742200% | 2.25052 |
| Density volume L1 | 0.0271766% | 0.0114011% | 2.38368 |

The128 thermal budget closes against the independently extracted snapshots: pressure work −1.070869823e−5, explicit compatible-time correction −2.681289361e−10, and net internal-energy increase +0.0001861419592. Positive accumulated heating still comes principally from the acoustic dissipation term; refining reduces it, but this does not prove a particular initiating error or exclude indirect time-discretization effects. It does not justify adding test-specific damping. The actual-cell momentum quadrature16→32 change is3.48e−9 absolute; reconstructed polygon areas agree with stored volume within9.67e−14 relative.

A finite CPU-only postprocessor now waits for the actual128 t=1.2 hydro/assessment completion, then applies the existing thermal-budget audit to all segments and saves the same64/128 comparison. It subsequently waits for glass Sedov completion and renders regular/glass density and contour comparisons with identical color/axis ranges and no smoothing. Live dependency checks and six-hour wait limits apply; failure never triggers a physics restart. Final glass figures require visual inspection.


## Current128 endpoint completed and compared at t=1.2

The current128 evolution and CPU polygon/thermal assessments have completed. The final segment accepted2401 steps with no retries; fixed masses, positive U/V and total area1 remain intact. Endpoint energy change is−4.40153e−15 relative, with maximum final-segment volume mismatch9.96667e−12. Source hashes and report paths are in `gresho_glass_evidence_20261109_v1/gresho_t12.json`.

| Metric at t=1.2 | Current64² | Current128² |
|---|---:|---:|
| Kinetic-energy loss relative to initial kinetic energy | 4.33362901 | 1.50133886 |
| Integrated acoustic heat | 0.00381391724 | 0.00132113339 |
| Cell-integrated momentum L1 (%) | 4.08223137 | 1.94912028 |
| Generator velocity L1 (%) | 5.00050077 | 2.39134772 |
| Density volume L1 (%) | 0.103476564 | 0.0392377751 |

Refinement reduces accumulated heat and kinetic-energy loss by a factor2.887 and momentum error by2.094. Spatial resolution and time cap were both refined; this is not a pure spatial-order measurement. The128 positive heat budget remains dominated by acoustic dissipation (+.001321133389), offset by pressure work(−.000063278378) and a very small explicit compatible-time correction(−1.54904e−9). This is a budget attribution, not proof of the initiating cause.

Compared with the historical128 method at the SAME t=1.2, current cell-momentum L1 is1.949120% versus1.890671%; current generator velocity L1 is2.391348% versus2.378210%. Density L1 is.0392378% versus.0396969%. Actual-polygon azimuthal-velocity panels were inspected and remain closely similar, including a slightly noncircular annulus/core. Thus common volume/work coupling is not demonstrated to improve Gresho accuracy by itself. Do not call this a pure inverse ablation because the historical code differs in other respects.
