# Paper priority and evidence register

The first astronomy-paper draft has priority. The second paper starts only
after that draft exists. Preparing for a fluid-mechanics journal means retaining
the following evidence during the current work; it does not authorize starting
a second manuscript or adding a separate campaign that delays the first draft.

## Evidence to retain now

| Question | Required record | Current limitation |
|---|---|---|
| What method produced each result? | Source and executable hashes, active options, IC/BC, resolution, target and achieved times, restart lineage | Successful results from different candidates cannot be pooled as one common method |
| Is the update physically consistent? | Fixed retained masses, momentum/energy including boundary transport and contact losses, geometric volume versus integrated volume, thermal-energy decomposition | Ledger closure alone does not establish accurate thermodynamics |
| Are the observed structures physical? | Actual unsmoothed cell polygons, pressure/velocity as well as density, shock and wall profiles, wake circulation and perturbation growth | A smoother wake is not automatically more accurate; seeded startup is not a growth measurement |
| Does the solution converge? | Separate spatial and temporal refinement at matched IC/BC and final time; independent reference uncertainty | Local convergence does not certify the whole domain or every boundary event |
| Is the method practical? | End-to-end geometry, weight, RHS and retry cost at matched accuracy and hardware | Helper speedups and shared-GPU timings do not prove the <=2x GIZMO objective |
| What Reynolds number is being studied? | Explicit physical viscosity and matched Re/Mach/BC when relevant | Current native/SU2 Euler comparison does not impose a finite physical Re; visual smoothness cannot determine it |

The working question for a later paper is whether matching geometric volume
evolution and pressure work reduces nonphysical thermodynamic errors in
fixed-mass Laguerre cells across shear, shocks and curved boundaries. This is
a research question, not a demonstrated universal result. New coefficients,
physical states or test-specific switches require a common derivation and
cross-test evidence. Earlier failed candidates must be searched before coding.

## Current result index and corrections

Let `D` denote `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1`
and `S` denote `D/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1`.

### New late-Cylinder evidence

`S/stage_volume_coordinate_v1/event_trace_v2` instruments the actual native
trajectory without changing its 13 extracted physical fields. At the retained
t=3.32844168159234 event, face 61601--62422 has almost equal reconstruction
increments around -0.01587. Their difference crosses zero, changing the old
ordered-jump factor from 0 to 1 and discontinuously restoring their common
normal-speed increment. This confirms the known limiter mechanism in this
specific late Cylinder event. It does not identify every prior error's cause.

`S/stage_volume_coordinate_v1/time400_v1/assessment_v1` records nine native
time-refinement runs from t=.01 to .0109765625 with 31,668 cells. Binary-exact
steps 2^-12, 2^-13 and 2^-14 give 4, 8 and 16 steps, without retries or radial
contacts. A preselected annulus containing 350 material IDs and nonuniform
fields shows observed maximum-norm orders about 2.00 for rho, P, U and velocity
with all three volume-update options. Whole-domain rho/V order remains about
1.00 near the wall. This is scoped temporal evidence, not a global second-order
or production-method qualification.

### Duplicated, rejected reconstruction withdrawn

The recent `S/continuous_normal_v2` free-trace formula repeats the historical
`D/continuous_jump` projection. The historical rejection was recorded in
`Docs/MATERIAL_RT_PORT_AUDIT_2026-11-06.md`, section "ordered-jump discontinuity
and common-speed counterexample". This history should have been checked first.
The repeated suite again fails Blast at t=.0004476789 of .038 and stalls Noh
at t=.005044567 of .6. Preserving an unbounded common reconstructed speed is
unsafe even when its relative normal jump is admissible.

The new late-Cylinder comparison completes its short interval, but cannot
override these shock failures. No common-method adoption follows. The opt-in
additions were removed from the five current C/CUDA/header files; 96 exact
repository-local source dependencies, hashes, executables and all results are
retained in `S/continuous_normal_v2/rejected_source` and its parent directory.
`withdrawal_verification/completed.json` confirms C/CUDA builds and byte-exact
Sod/Blast initial and final outputs against the previous default control.

The historical `D/bounded_jump` candidate jointly limits common and relative
speed by one projection onto stencil extrema intersected with the ordered-jump
strip. It completed all nine 1D tests and the recorded KH/shear/Sedov/Gresho
tests, but did not improve accuracy uniformly: Sedov64 rho/P L1 became
.11196878/.08705445 versus .10891516/.08447629. Reuse its mathematical helper
if appropriate; do not transplant the obsolete compact-velocity ABI or its
stationary-wall-only driver into the current Cylinder implementation. Current
prescribed exterior traces and moving image-wall work need explicit treatment.

### Additional swept-CFL continuation is terminal

`S/common_swept_cfl_v1/continuation_t5_v1` exited 99. Last accepted time is
2.8378217039438036; the most recent valid saved checkpoint is
`cylmat.007697.dat` at 2.7915000000000014. The next retried predictor reached
the "invalid generator at cylinder" check. The endpoint worker correctly
recorded failure and launched no restart or false t=5 analysis. The older
native800 baseline already completed t=5 and remains intact.

Subsequent logging-only replay reproduced all 133 accepted records and
identified the actual failure: the original reflection executes, but its
coordinate update rounds back to an inside point. Private V5/V6 variants
qualify the same reflection with a stored-coordinate predicate; the successive
entry/postcondition failures and current terminal outcome are recorded in
`Docs/CYLINDER_CONTACT_REPRESENTATION_2026-11-08.md`. These are numerical
boundary-representation findings, not physical wake instability or full t5
qualification. Preserve the original velocity law and contact-energy ledger.

### Positive-gap follow-on (latest)

The private common radial metric now preserves positive represented wall
clearance even when hypot rounds to R. C/CUDA geometry and actual native
failure/follow-on tests pass through2.839. Its later continuation toward5
was intentionally stopped after a local thermodynamic failure was measured,
as recorded below. All failed runs and the old t5 baseline remain preserved.
No complete-t5, morphology or canonical-adoption claim is made. See Docs/CYLINDER_CONTACT_REPRESENTATION_2026-11-08.md and its linked
positive_gap_v1 evidence and live/terminal process records.

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

## Immediate order

1. Diagnose the actual Cylinder boundary rejection and qualify its geometry,
   thermodynamics and morphology using the existing code and retained states.
2. Reconcile the common reconstruction with existing successful and rejected
   candidates; avoid repeating old experiments or accepting a local improvement
   that fails shocks elsewhere.
3. Revisit the remaining 2D comparisons, then qualify full 3D fluid evolution.
4. Use demonstrated results and explicit limitations in the first draft.
   Begin the second manuscript only after that draft is available.

### Local entropy/volume follow-up

Read-only same-trajectory diagnostics and actual smaller-step tests now
separate local entropy loss from global energy closure. The selected cell's
RHS heat is positive, but rebuilt volume differs from integrated Vdot and
dominates its negative entropy change. Decreasing the caps below the actual
event step reduces peak volume discrepancies1.4284% ->0.46284% ->0.24851%.
The finest short interval has positive net entropy change, although individual
steps remain nonmonotone and no practical correction is adopted. See
Docs/CYLINDER_LOCAL_ENTROPY_2026-11-08.md and positive_gap_v1/entropy_trace_v1.
This is a local mechanism diagnosis, not full Cylinder or2D qualification.

At the exact maximum-error endpoint, the existing first-order derivative
applied to all neighbouring contact position changes accounts for99.5158%
of the selected volume residual. The selected cell has no direct contact.
This is scoped sensitivity evidence, not a finite-jump proof or final method
adoption. See the local entropy record for replay identity and limitations.

### Finite-precision stopping evidence

The existing volume inverse, reconnected to current native contact geometry,
passes disabled byte-exact parity but both enabled modes stop before any
accepted step. New measurements show that rounding W+=delta predicts the next
near-contact volume residual: 3.476911027e-9 predicted versus 3.476910994e-9
observed, whereas the requested linear correction residual is about 1e-24.
This distinguishes a storage/convergence issue from a new physical instability.
It does not qualify the later contact event or the full Cylinder solution.
The detailed record includes failed-run costs, 100-digit independent checks,
the invalidated initial inline diagnostic and authoritative stored-W analysis.
See Docs/CYLINDER_LOCAL_ENTROPY_2026-11-08.md, final section.

Preparing for the later fluid-mechanics journal therefore includes retaining
negative results and separating physical-model error, spatial/time truncation,
boundary-condition effects, and finite-precision stopping errors. Work on the
first astronomy draft and its existing validation remains the priority; the
second manuscript begins only after that draft exists.

### Completed volume/pressure-work comparison

The existing volume inverse with a separately identified finite-precision stop
now completes the same native interval through2.839. The selected cell's
entropy proxy changes +2.2356% (endpoint) or+2.2335% (stages), versus-1.9515%
without the inverse; geometric volume discrepancies stay below3.48e-9. This is
evidence for the later paper's physical question, not whole-method validation.
The inherited earlier error and stationary-wall accuracy remain unresolved.
See Docs/CYLINDER_CAPACITY_PRECISION_2026-11-08.md for dimensional scope, speed
limitations, reproduction and the private symbolic-analysis reuse experiment.
The first astronomy draft retains priority; no second manuscript is started.

The private native solver now reuses unchanged CSR structural analysis while
retaining full nonsymmetric factorization and raw-geometry residual checks.
Both mode comparisons complete with density/pressure differences below7.4e-10
against fresh solves. Repeated analysis cost decreases materially; the shared
GPU totals do not yet establish the GIZMO performance requirement. The first
astronomy draft still takes priority over beginning a second manuscript.

### Independent earlier contact-volume sign check

The same native capacity inverse now also completes the saved t=.975->.9765
contact event (ID64062), whose previous geometric entropy defect was positive.
All-cell volume residual falls from1.6833% to<=1.48e-13; the selected entropy-
proxy increase falls from0.6893% to0.00999--0.01075%. Retained masses are exact,
and all three runs exit0. Stationary-wall physical energy errors remain.
This supports the common mechanism across two opposite-sign local events;
it does not qualify whole-domain physical entropy, the inherited ID63222
history, final morphology or the full method. First-draft validation remains
the purpose; no second manuscript or separate campaign is started. Evidence:
Docs/CYLINDER_CAPACITY_PRECISION_2026-11-08.md, final section.

### Before-loss restart and matched-reference comparison

The same saved t=2.5415 state was continued with the existing volume inverse
off/endpoint/stages to t=2.65. The known cell's density error against an
independent same-time fine reference falls from178.7% to2.43--2.47%; its
large entropy-proxy loss is prevented. Both reference mesh resolutions and
actual native cell-average projection are retained. Whole-cohort improvement
is modest compared with reference mesh sensitivity. Stationary-wall energy,
boundary shape and a downstream cell-average entropy decrease remain open.
The next controlled interval reuses these t=2.65 checkpoints, without a new
IC. See `Docs/CYLINDER_LOSS_PREVENTION_2026-11-08.md` and its raw evidence.
This continues first-draft validation while retaining evidence for the later
fluid-mechanics paper; it does not start the second manuscript.

## Physical justification and necessity (2026-11-09)

Every proposed change must answer both questions: is it physically or
mathematically justified, and is it necessary in the simplest common method?
Improving one benchmark is insufficient. Prefer removal/replacement of an
inconsistent operation to adding a test-specific sensor, coefficient or
physical state. Preserve the successful Kepler baseline as evidence; do not
silently change it or count different candidates' successes as one method.

The current volume inverse is a numerical compatibility candidate: the same
discrete pressure-work volume target is realized by changing W, without
repairing U, momentum, positions or cell masses. Its actual native support is
single-rank 2D Cylinder. The anchored full-J implementation cannot be used
unchanged in a periodic box: constant-weight gauge freedom and compatibility
must be handled, then the original full-row residual checked. Existing
periodic gauge/direct-solve and static Gresho inverse evidence can be reused;
they are not stage-by-stage general2D evolution qualification. 1D moving
faces, general2D, MPI and3D are not validated by this Cylinder improvement.
See Docs/COMMON_METHOD_NECESSITY_2026-11-09.md for evidence and the rejected
offline structural-envelope performance experiment.
