# Cylinder: preventing the inherited local density/entropy error

The first astronomy draft and its existing validation remain the priority.
This record preserves evidence for a later fluid-mechanics paper; it does not
start that manuscript or a separate simulation campaign.

## Controlled comparison

Let `S` be `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1`;
`A=S/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1` and
`F=A/prevention_t265_v1`.

All three runs restart the SAME saved state at t=2.5415000000000014, before
the known large local entropy loss, and complete t=2.65. The initial
800x160 grid has 126,696 current owned cells. The executable is identical:
SHA256 `7d61cb4a6a51b0252a92db80f76c13bdc4b9f9b54c1a4d6a8ba589d14a891ccf`.
The only physical-algorithm difference is the existing volume inverse:
disabled, applied at the endpoint, or applied at stages and the endpoint.
Source/checkpoint/companion hashes and full environments are in each run's
`lineage.json`, `environment.json` and the three `*_plan.json` files.

The inverse solves for Laguerre weights to realize the integrated cell
volumes. It changes neither stored mass nor internal energy by a repair.
Original radial reflection, inward-normal velocity removal and the explicit
contact kinetic-energy ledger are retained. This is an Euler/slip-cylinder
test, not a finite physical Reynolds-number comparison.

## Local result, with independent references

The previously identified material ID63222 is followed without reselecting a
favourable cell. A=P/rho^gamma uses the effective stored gamma. A/start ratios
refer to this restart, not to the IC; the initial state already has A/A_IC=1.641.

| Quantity | Inverse off | Endpoint | Stages + endpoint |
|---|---:|---:|---:|
| Accepted steps | 279 | 244 | 242 |
| A_final/A_start, ID63222 | 0.419499 | 1.006466 | 1.002432 |
| Final density, ID63222 | 10.989997 | 4.017425 | 4.015130 |
| Final pressure, ID63222 | 1.578800 | 0.925818 | 0.921370 |
| Maximum all-cell step volume residual | 0.244507 | 3.202e-11 | 6.196e-10 |
| Density error versus fine cell-average reference | +178.69% | +2.47% | +2.43% |
| Pressure error versus fine cell-average reference | +67.87% | -0.46% | -1.01% |
| Density error versus peer cell-average reference | +193.77% | +7.54% | +7.69% |
| Pressure error versus peer cell-average reference | +70.80% | +2.84% | +2.37% |

In the off run, the selected density transiently reaches 30.81884. The full
unsmoothed history, including this maximum, is shown in `F/thermodynamics.png`
and `.pdf`; no visual clipping or smoothing is used. Both enabled modes meet
the predeclared 1e-8 per-step and 1e-6 cumulative volume budgets. The largest
cumulative sum of all-cell maxima is 1.878e-9. Separately certified
finite-precision stops remain distinguished from normal strict convergence.

Independent SU2 h=.025 and .05 fields were continued from their existing
BDF2 checkpoint pairs to the SAME t=2.65, without a new IC. Both complete
without nonphysical recovery or inner convergence failure. Native fields are
compared with reference integrals over the ACTUAL native polygons, rather
than point samples. Maximum projection coverage error is below 2.75e-14.
The references are in `cylinder_reference_su2_v1/matched_h0025_t265_v1` and
`matched_h005_t265_v1`, alongside the existing Cylinder comparison directory.

This provides evidence that matching actual and integrated volume prevents
this large local thermodynamic/density error. It does not establish universal
entropy preservation or qualification of the entire method.

## Whole-region and boundary limitations

The prespecified upstream cohort contains 1,834 material IDs selected using
the common initial position (x<2.5, distance from (2.5,2)<.7), not final error.
Density relative L1 changes 2.3650% -> 2.2824%/2.2825%; pressure relative L1
changes .9739% -> .8979%/.8985%. Fine/peer reference sensitivity is about
1.782% for density and 1.641% for pressure. Those modest aggregate changes
are smaller than the reference mesh difference; do not claim a resolved
whole-region accuracy advantage from them.

All 126,696 final polygons in every mode are convex, with no cylinder
intrusion, and reconstructed volumes match the saved values. Nevertheless,
the excluded obstacle area remains about .7917799 instead of pi/4=.7853982;
the largest boundary radius is about .526815 rather than .5. Volume consistency
has not removed this boundary shape error. The actual wall-pressure spike
falls from 2.12754 to 1.02528, but pressure dips remain. See the actual native
face-pressure/reference comparison in `F/wall_pressure.png` and `.pdf`.

All retained material masses are exact. The entropy population is the
126,056 IDs shared by the initial and final snapshots, not all 126,696
snapshot IDs in this open-boundary problem. Its lowest delta log A changes
from -.868694 at ID63222 to -.008524/-.008528 at ID60066. ID60066 is downstream
(final x~3.34737, y~1.89812), outside the upstream reference cohort and not a
wall-face owner. Its detailed RHS entropy budget was not recorded. Negative
cell-average entropy-proxy change alone does not prove a physical violation;
prior exact-boundary-flux diagnostics already demonstrated this distinction.

Stationary-wall physical energy defects remain -2.39369e-5, -2.25409e-5 and
-2.39256e-5, respectively. These combine image-wall work and original contact
loss; the ledger residual itself is below 1e-14. Ledger closure must not be
reported as eliminating stationary-wall work or contact losses.

Elapsed times are 436.8/687.6/1007.8 seconds on a shared Ada GPU. Concurrency
differs, so these are operational costs, not matched-accuracy performance
measurements. The <=2x GIZMO objective remains unverified.

## Reproduction and continuation

`F/completed.json`, `F/assessment/`, `F/*_geometry/`, `F/*_reference/`,
and archived run/build/analysis scripts preserve the result. The observer
calls current C/CUDA RHS without modifying it and exits before any evolved
step; immutable-state checks and exact rebuilt-volume agreement pass.

The three calculations now continue from their own verified t=2.65 snapshots
to the existing historical bracket endpoint t=2.7915000000000014, under
`F/continuation_t27915_v1`. Only the target/output schedule changes; source
physics, executable and options are inherited. Existing outputs are retained.
Live versus terminal state must be read from process records, not this
description. No second manuscript, 3D evolution or canonical-method adoption
is implied by this controlled continuation.

Historical review before further changes found that U-only adiabatic splicing
leaves unmatched weight work; weighted KDK uses an auxiliary kinetic-energy
state and is restricted to smooth periodic tests; frozen acoustic-face flow
fails positivity at some large steps. None is a proven drop-in Cylinder fix.
Do not repeat these candidates under new names or change time integration
merely because a cell-average entropy proxy decreases.
