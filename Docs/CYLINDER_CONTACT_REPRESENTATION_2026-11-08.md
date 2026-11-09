# Cylinder contact representation audit (2026-11-08)

The first astronomy draft retains priority. This is verification evidence for
that work and possible later fluid-mechanics analysis, not a second manuscript.
The older native800 t=5 baseline remains preserved separately.

## Confirmed original failure

The logging-only replay from the retained t=2.7915 state reproduced all 133
accepted records exactly before the additional swept-CFL branch failed.
At t=2.8378217859584858, predictor phase2, ID63220 has represented squared
wall clearance -5.362203420604374e-18. The original radial reflection executed
and removed inward velocity, but its position update rounded back to the same
inside point. sqrt is just below R and hypot equals R. Periodic wrapping did
not change this particular failing point. Source and rational proof are in
`Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/terminal_audit.json`.

## Narrow private candidate and actual qualification

The same original reflection is evaluated as
`x_reflected = x - 2 D d/[r(R+r)]`, with `D=d.d-R^2<0` and `d=x-c`.
Compensated arithmetic retains the penetration and each output coordinate is
rounded once outward to enclose that reflected value. No physical clearance,
repeated displacement, new velocity law or internal-energy repair is added.

The first native guard failed because relaxed C arithmetic tested
`n*r_out` after optimizing `(center+n*r_out)-center`, rather than inspecting
the stored coordinate. The compiled actual contact fixture proves this.
V5 separates the stored-coordinate postcondition into a precise, non-inlined
function. Its actual fixture preserves vx,vy,U and contact delta KE exactly.
Two actual calls then match independent directed high-precision references,
but the run next fails at phase3 ID63227: a genuinely inside point has
sqrt==R and bypasses the old collision-entry condition.

V6 uses the same compensated squared-gap sign to enter and check reflection.
Both actual failure inputs pass compiled contact fixtures; the original first
case is unchanged from V5 and the second uses the original inward-velocity
removal law once. C/nvcc-host/CUDA tests have 3,192 identical helper outputs,
3,184 accepted and 8 conservative rejections. All accepted outputs match the
independent directed reference and have positive exact-rational clearance.
The unchanged normal one-step test retains all 13 extracted fields exactly.
These are limited numerical qualifications, not an all-input predicate proof.

## Latest native outcome

V6 is terminal (exit 99) at t=2.8378227218785139, phase 2, ID63230. The represented failing point is outside; exact squared gap=1.7998348912037156e-17. It recorded 140 accepted steps, preserved 133/133 original accepted records exactly, and invoked the representation fallback 6 times. New checkpoint filenames: ['cylmat.007823.dat']. Check restart metadata and validity before reuse.

The existing Astep=.04632 output mechanism is used to save the first accepted
state after the pre-event threshold. Astep is not part of the dt ceiling;
HYDRO_TSTOP is unchanged. This avoids cutting the step schedule merely to obtain
a diagnostic restart. Source hashes, commands, traces and immutable binaries
are in `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/rounding_candidate_v6`. Structured candidate outcomes are in
`Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/contact_candidates_summary.json`.

## Remaining distinctions

Exactly zero represented gap is not solved by this reflection helper. A true
positive gap whose hypot rounds to R is also unsupported by the current image
chain: radial_split_point computes displacement from rounded R-r, making the
image coincide; the self-face-speed and strict padding checks reject it.
High+low coordinates preserve a nonzero image displacement only after that
displacement has been computed. They cannot recover a gap already lost in r.
Removing those checks alone creates degenerate geometry.

The 8 synthetic rejections remain documented; a rejection currently terminates
rather than silently moving a point farther. Keep true contact, positive-gap
representation, and actual penetration distinct. No whole-method adoption,
global-order, MPI, morphology, full3D-fluid or matched-GIZMO-speed claim follows.

## Verified pre-event restart and next implementation boundary

The new checkpoint `rounding_candidate_v6/event/cylmat.007823.dat` is verified
at t=2.8378200984135984, 800x160, 126,696 cells. Its reader, unique IDs, finite
fields, positive U/V and all three file hashes pass. The 126,376 retained IDs
have exactly unchanged mass; 320 entered and 320 left through the open boundary.
Do not mistake this documented ownership change for variable material mass.
The minimum rounded radius is .5000000000000033. The accepted record is the
126th record of the unchanged prefix. See checkpoint_audit/completed.json.

The next question is strictly positive represented clearance whose rounded
hypot is R (actual ID63230), not another penetration. A candidate must keep the
fluid coordinates unchanged and carry `(d.d-R^2)/(r+R)` into image displacement
before it is lost. Its sign, high+low geometry, reflected vectors and face
speeds need a common C/CUDA check; simply accepting r==R is insufficient.
Exact D=0 remains a distinct geometric degeneracy. No such image-chain change
has been implemented or qualified yet. Restart from the verified near-event
state; do not replay from t=2.7915 or initialization without a specific reason.

## Positive-gap implementation and actual continuation

The next-action paragraph above is superseded by the implemented private
`positive_gap_v1` common radial metric. It preserves the positive gap in
image construction and the same side in self-face speed, without moving an
exterior fluid point. Its C/nvcc-host/CUDA high-precision geometry checks and
existing support certificate checks pass. The actual failure interval and
follow-on both complete, reaching t2.839 with fixed masses and positive U/V.
The additional swept-CFL candidate is now continued from2.839 toward5; it is
not yet a t5 result or a canonical adoption. See
`Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/
cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/
production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/README.md`
for hashes, numerical evidence and the live/terminal state source.

### Thermodynamic follow-up and intentional stop (latest)

The positive-gap continuation was intentionally terminated after 557 accepted
steps, last accepted t=2.8708767579956156. Its native exit is -15 (SIGTERM),
not a new physical abort or target completion. The verified saved restart is
t=2.839; later accepted log states are not saved restart states. All outputs,
failed experiments and the older separately completed native800 t5 baseline
are preserved. This stops only this additional trajectory; root-cause work
continues.

For retained material ID63222, A=P/rho^1.4 decreases from the actual IC value
0.08000000559 to0.02803433007 at t=2.8378200984 (35.0429% of the IC value).
This error already exists before the new positive-gap geometry change. At
t=2.839 it is0.02681755223 in the coarse follow-on. Global energy-ledger
closure does not qualify local thermodynamic evolution.

Same-source adaptive-cap comparisons to2.839 give maximum local volume-step
quadrature residuals4.7170%,4.4034%,4.4034% for caps.0005,.000125,.0000625.
The latter two maxima occur at the same actual dt=2.513864e-5, below both
caps. These data therefore do not establish irreducible error or temporal
convergence. The coarse/finest ID63222 endpoint density differs5.888%, while
whole-domain material-mass-weighted relativeL1 differs only7.85e-6.
Next diagnostics resolve local H=Udot+P Vdot, actual finite-step U/V changes,
and caps below that event's actual dt, preserving the same source checkpoint.
No thermal-energy repair, new physical parameter or method adoption follows.

Evidence in positive_gap_v1: continuation_t5_intentional_stop.json,
continuation_t5/process_result.json, entropy_audit/completed.json,
time_refinement/assessment/completed.json. The first astronomy draft remains
the priority. These retained diagnostics prepare the later fluid-mechanics
paper without beginning another manuscript or unrelated test campaign.
