# Rotationally symmetric Sedov source control — 2026-11-09

## Reason and scope

The64²/128² square-source orientation controls changed the fourth-mode direction when the source rotated45°. That implicates the source and its projection, while leaving grid coupling and operator errors unresolved. The next controlled change is source geometry, with the same common C/CUDA hydro binary, fixed masses, weights, lattice, EOS and boundaries. This is first-paper2D validation, not a3D campaign or a new viscosity/force/weight correction.

The user asked whether256³ would reduce anisotropy. That is a3D resolution; the current2D comparison is256². Resolution can reduce discrete projection errors but cannot remove anisotropy of a fixed physical square. The aligned source's fourth contour amplitude increased64²→128², so larger N alone is not established as a cure.

## IC implementation and dimensional comparison

Existing `Exam/NativeU64/unified_initial2d.c` now optionally accepts `GFS_SEDOV_SOURCE_SHAPE=disk`. A physical `GFS_SEDOV_SOURCE_WIDTH` must be explicit. The disk radius is width/sqrt(3), so both the disk and square have the same continuous energy-weighted second radial moment width²/6, as well as the same total excess energy1 and center(.5,.5). This matches a declared spatial extent measure, not the full radial profile or support area. At width1/16, radius=.03608439182435.

Analytic rectangle/disk intersection integrates the source into each initial rectangular cell, using a long-double quadrant primitive. No angular polygon approximation, numerical smoothing, renormalization, generator shift or hydro change is used. Invalid shape, missing physical width and nonzero effective disk orientation are rejected. Default and explicit-square paths retain original arithmetic. Running native binaries remain unchanged.

## Completed verification

Before/after IC writers reproduce default KH, Sedov, Gresho, periodic shear, periodic Sod and Noh checkpoint files bitwise. Explicit square0°/45° reproduces archived source-control checkpoints bitwise.1000 off-grid random rectangles agree with independent adaptive vertical-chord integration to max area error5.56e−17. Disk32²/64²/128²/256² deposited cell energies agree with independent quadrature to1.57e−17 absolute, total excess energy differs from1 by≤2.23e−16, and mass/volume/zero velocity/zero initial weight are unchanged. Fine-to-coarse deposited thermal energy differences are≤2.78e−17; mass and volume coarsen exactly. No evolution at256² is claimed.

A circular continuous source still projects onto square cells. Its projected H4/R4 moment is+.0420762 at64²,+.0134844 at128²,−.00143642 at256²; the continuous value is0. These moments quantify initial spatial projection, not the final contour Fourier amplitude. Projected second radial moments also approach the matched continuous value; projection error must remain part of the interpretation.

## Execution and output contract

Detached finite supervisor443022 runs the archived `run.py`. At launch it is waiting for successful current Gresho completion and actual exit of its final native process. This reuses the third Ada GPU slot and never stops/restarts cylinder or RT.64² then128² evolve to t=.1, C=.4, dt caps.000125/.0000625, unchanged portable binary SHA11a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e. Fixed mass/positive state/area/energy and endpoint-volume≤1e−11 gates remain.

Effective parameters are written to `run64/params.dat` and `run128/params.dat`. Output interval.05 retains a midpoint and final state, at most4 checkpoints including input per case, total budget57,017,248bytes; launch free space149,288,197,816,320bytes. CPUs12,13 evolve;20,21 analyze saved power-cell geometry. Existing Sedov analysis measures reference discrepancies and complex contour modes. Original square results remain the controls. Point-source analytic curves are not exact finite-disk solutions; endpoint completion is not a morphology pass or a GIZMO-speed certificate.

Evidence: `source_change.json`, `initial_audit.json`, `plan.json`, `launcher.json`, `status.json`, frozen before/after writers and run scripts in `common_sedov_disk_20261109_v1`. No commit/push, no3D evolution, no manuscript revision.


## Completed64 evolution and three-source comparison

The disk64² case completed t=.1 with814 accepted steps,14 retries,147.492s shared-device solver wall time. Exact masses, positive U/V and total area1 remained; max volume mismatch9.83062e−12 and relative global energy change+3.84015e−13. Actual128 native process449550 is running; no128 disk endpoint is yet claimed.

All three64² rho=2 contours have256/256 rays. The aligned-square/disk comparison reduces fourth-mode amplitude1.11710%→.135684% (about88% reduction of this mode), but total angular RMS only1.19662%→.903961%. The fourth mode accounts for approximately43.6% of the aligned contour variance and1.13% of the disk contour variance. This does not mean88% of all morphological error has disappeared. Other angular components, cellwise contour sampling and interior structure remain.

Ideal point-source density L1 mismatch is10.7898% aligned square,12.0830% rotated square,10.9362% disk; pressure mismatch8.22777%,11.4706%,9.57764%. A more circular contour alone therefore does not establish lower radial error. The source profiles differ from a point explosion and from each other even though energy and second radial moment are matched.

`figures64/density.png` and `angular.png` compare all three without smoothing or selective omission, using actual power cells and a shared color range0–5. Visual inspection confirms the disk result retains localized density deficits and diagonal structure behind the front. The source-geometry contribution is reduced in the predeclared fourth mode; full Sedov morphology is still unresolved. Wait for128 before interpreting a resolution trend; no additional force/viscosity/limiter patch is justified by this64 result alone.


## Completed128 disk and resolution comparison

The128² disk evolution reached t=.1 in1662 accepted steps,26 retries and643.392 shared-Ada solver seconds. Exact fixed masses, positive U/V and area1 remain; max endpoint volume mismatch9.9442e−12, relative energy drift−2.2325e−14. Saved-state polygon analysis and three-source plots completed; `figures128/density.png` was visually inspected. No smoothing or selective omission is applied.

For the disk, rho=2 contour angular RMS decreases64²→128² from.903961% to.398678%, while fourth-mode amplitude changes.135684%→.129706%. Density point-reference L1 decreases10.9362%→6.39036%, pressure9.57764%→4.28370%. At128 the aligned-square contour RMS is1.17417%, so the isotropic source reduces this diagnostic substantially. These are threshold/profile diagnostics, not a proof of the exact finite-source solution or formal convergence order from two resolutions. Fine radial streaks and local density structures remain visible behind the front. Hence the user's glass-layout proposal addresses a remaining distinct question; source shape and layout must be controlled separately.

The disk128 result and analysis released the third Ada slot to current Gresho128 as scheduled. No regular256² hydro or any3D evolution has been launched.256² disk IC exists only as a verified initial-deposition artifact.
