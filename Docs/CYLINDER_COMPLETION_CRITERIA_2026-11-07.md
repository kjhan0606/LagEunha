# Cylinder completion criteria

User priority: complete cylinder first, then revisit other 2D tests. Completion
means physical results at least comparable to established work. A short stable
run, a clean energy ledger, or a rendered density plot is insufficient.

## Physical target and provenance

Retain the manuscript target as the existing project baseline: inviscid
compressible flow, rho=1,U=1,P=.08,gamma=1.4 (Mach2.988...), cylinder R=.5 at
(2.5,2), domain20x4,400x80 initial sites, endpoint at least t=5. Verify the
actual historical input and boundary conditions before treating the prose as
an executable specification. Do not relabel Mach.1 native boundary diagnostics
as this benchmark. Preserve old calculations and restart when compatible.

Locate an independent published/reference calculation with matched Mach,
gamma, obstacle/domain, inlet/outlet and transverse boundary conditions.
A channel with slip side walls, periodic sides and an unbounded exterior are
not interchangeable. Reynolds-number-dependent viscous shedding cannot be
used as the reference for an inviscid calculation. Compare matched
nondimensional time Ut/D and h/D, not raw cell counts alone.

## Required evidence

- Stationary obstacle: measured normal wall motion and wall power; penetration,
  geometry and boundary errors under resolution refinement. No position fixes
  masquerading as a successful no-steering material update.
- Open boundaries: consistent exterior geometry and state at both SSPRK2
  stages, explicit incoming/outgoing mass,momentum,energy ledger; owned masses
  unchanged between documented ownership-entry/exit events.
- At least the target time and spatial resolution above; representative bow
  shock and developed wake at matched times, not only the initialization phase.
- Density/pressure morphology: bow-shock stand-off, full shock contour,
  surface pressure and downstream cross sections compared with an independent
  matched reference. Also measure symmetric upstream structure without forcing
  physical unstable wakes to remain symmetric.
- Resolution/time-step refinement establishes error, rather than merely
  changing the colour scale or smoothing images. Quantitative acceptance bands
  must use reference uncertainty and matched reference resolution; no post-hoc
  thresholds chosen to bless this implementation.
- Wake vortices can be physical Kelvin--Helmholtz roll-up of the shear layers.
  A smoother reference image is not automatically more accurate. Compare
  velocity/vorticity with the density structures and check shear-layer growth
  under matched spatial/time refinement. Pointwise density differences alone
  do not distinguish physical instability, differing vortex phase, numerical
  diffusion, and spurious cell-scale fluctuations. Preserve the distinction
  between the present transverse-periodic calculation and published wall-
  bounded Mach3 cylinder flow; do not impose a viscous Karman-street template.
- The executed SU2 comparison uses Euler/HLLC/MUSCL/Barth--Jespersen, with no
  imposed finite physical Reynolds number. Stronger scale-dependent numerical
  dissipation could make its effective Re lower and suppress smaller wake
  structures, but this is not quantified by visual smoothness. The same-h=.05
  dt-halving comparison changes wake densityL1 by2.586e-5 and velocityL1 by
  6.218e-6; that does not isolate spatial dissipation. Do not label SU2 as a
  low-physical-Re solution or calibrate the native solver to its smoothness.
  Executed settings and limits: cylinder_reference_su2_v1/
  reynolds_interpretation_v1/audit.json.
- Positive thermal energy/cell volume, convex valid cells, fixed material
  masses and conservation including boundary transport. These checks remain
  necessary but cannot replace morphology.
- Runtime including weight/geometry work, and the thread's <=2x GIZMO objective
  at matched accuracy where the same physical benchmark is implemented.

## Status

Latest endpoint-volume update (2026-11-09): the existing private candidate
now completes to3.00025 from the before-loss2.5415 restart lineage. All126,696
polygons are convex with no intrusion; fixed retained masses and positive U/V
pass. Same-time SU2 pressure drag differs0.0875%, but wake cell-average
density/pressure differences3.049%/4.467%, approximate wall shape and nonzero
wall work remain open. See Docs/CYLINDER_MATCHED_T3_2026-11-09.md for complete
scope and evidence. A same-binary continuation from3.00025 toward3.50025 is
the next validation interval; this is not a newly qualified t5 or a common
cross-dimensional method. No second manuscript has been started.

Historical contact-representation update (2026-11-08): the earlier positive-gap
trajectory stopped at accepted2.870876758 after a thermodynamic error;
its saved restart2.839 and earlier failures are preserved. The new lineage
above starts before the known loss. These are different trajectories.

Retained earlier baseline (2026-11-08): the retained native800 baseline completed t=5.
All126,696 owned cells pass the endpoint polygon/nonpenetration check;
retained masses are exact and U/V remain positive. Whole-lineage physical
relative energy defect is -1.72988e-5, with ledger closure -6.70e-14.
Complete-near native400-to800 comparisons against the same finite SU2
reference improve interior density difference3.31% to2.22%, wake density
10.47% to7.87%, and pressure drag difference1.94% to.61%. Different late dt
caps prevent a pure spatial-order claim. Wall-area excess and local nose
pressure fluctuations remain. The baseline does not include the new swept CFL.

The complete11,793-step history contains30 steps with >1% maximum local
volume-increment quadrature residual and a44.2627% maximum at t3.328923928.
The1% level is descriptive, not a new acceptance threshold. Separate local
swept-CFL tests reduce the known extreme events, but have not qualified the
whole t5 trajectory. Small global energy error cannot excuse a local defect.

Both native800 and SU2 saved velocities exhibit opposite-sign upper/lower
wake circulation. This supports the shear-layer interpretation, not a KH
growth-rate measurement or certification of every density feature. Physical
KH must not be damped merely to resemble the smoother reference. The common
swept-CFL pilot from t=.27025 to.29025 completed: both settings retain exactly
the same13 extracted fields,41 steps, with an observed shared-GPU cost ratio
1.078. The bound is nonlimiting there; this is routing/cost evidence, not
improved accuracy. The additional adaptive branch in
common_swept_cfl_v1/continuation_t5_v1 is now terminal (exit99): last accepted
t=2.8378217039438036, latest valid saved checkpoint t=2.7915000000000014.
A retried predictor failed the "invalid generator at cylinder" check.
A subsequent logging-only replay preserves all 133 accepted records and
records the actual inside point. The original reflection executes but loses
its coordinate update on storage; a second event also exposes a rounded entry
test. Private representation candidates and their precise scope are recorded in
`Docs/CYLINDER_CONTACT_REPRESENTATION_2026-11-08.md`. They do not relax the wall
gate or certify the complete swept-CFL trajectory. The original endpoint worker
recorded failure and made no automatic restart. This additional-run failure does not replace the retained
baseline t5 result above. Earlier errors in the inherited state are not erased.
No3D fluid evolution is authorized by
these endpoint results; full2D morphology, MPI and matched-accuracy cost remain
open. See production_v1/spatial800_v1/handoff_t32375_v1/automatic_extent_v1/
endpoint_t5_v1/README.md and assessment.json for exact evidence.

Source-only speed update (2026-11-08): the common geometric-CFL helper now
reuses the same robust distance rather than calculating it repeatedly.
All1,209 C/CUDA checks and the short KH/RT/Sedov/Noh plus early/late Cylinder
states and time bounds remain exact. Isolated helper cost decreases about31%;
this is not a whole-application or matched-GIZMO claim. The live continuous
candidate keeps its earlier immutable executable. Qualification and hashes:
production_v1/spatial800_v1/common_swept_norm_reuse_v1/adoption.json.

Static endpoint-volume result (2026-11-08): at the retained t3.328923928 event,
the actual integrated SSPRK2 target volumes are all positive. A private frozen
native-geometry inverse, using the full image/alias/exterior volume Jacobian,
reduces the44.2627% residual normalized by initial volume to1.695e-12 in four
Newton corrections. Positions, velocities, masses and energies remain exact;
all126,696 polygons retain positive area/convexity and pass wall nonpenetration.
This proves a static weight solution exists for this endpoint, not a valid or
fast projected time integrator. The10.12s diagnostic cost is not a production
or matched-GIZMO measurement. No change was applied to the live swept-CFL run.
See spatial800_v1/endpoint_volume_target_v1/README.md for qualification/limits.

### Earlier observations (superseded by the endpoint above)

Earlier update (2026-11-08): the t3.2375 inlet failure is reproduced and
isolated to the finite computational search bounds. Automatic expansion
of measured pure-halo shortfalls is qualified against fixed wider bounds,
physical padding controls and unchanged short KH/RT/Sedov/Noh states.
The 800x160 branch has resumed from t3.2385 toward t5; no t5 completion or
morphology qualification is claimed yet. See
`spatial800_v1/handoff_t32375_v1/README.md` under the cylinder production
research directory. Earlier status paragraphs below are historical.


Current native800 observation (2026-11-08): main target5 was not reached.
The last complete Cylinder acceptance is3.2373428361583056. The following
SSP step reached3.2375, but the inlet ownership handoff failed the finite
GPU search-region certificate for160cells (exit99). No t5 morphology has
been produced for this branch. A same-binary replay from the saved3.021
checkpoint is recovering a pre-event state; see
production_v1/spatial800_v1/handoff_t32375_v1. This is a boundary/search
failure, not a rejection of physical KH wake roll-up.

The isolated full-face swept-volume timestep candidate reduces a different
wake event's local volume residual5.14%→.158%, with final selected density
and pressure within.011% of a finer fixed-step calculation. It is still an
experimental build and is not in the failed main executable. Full-t5 time
convergence, wake morphology and matched-accuracy runtime remain pending.
See spatial800_v1/wake_volume_t2672_v1/README.md. Earlier observations below
are historical and do not imply that native800 is currently running.

Implementation update (2026-11-08): the redundant second full face-normal
preparation now refreshes walls only after the recovered wall speeds change.
Frozen C/CUDA and actual native short evolutions remain byte-identical, with
an additional canonical unprofiled build/contact check. Source provenance is
in production_v1/wall_normal_refresh_v1/adoption.json. This reduces duplicate
work without changing the method or satisfying any outstanding physics gate.
The ongoing native800 calculation keeps its earlier immutable executable;
it was verified live at t=1.651, without restart. This time is an observation,
not an endpoint. Short shared-GPU timing is not the matched-GIZMO speed test.

Latest refinement update (2026-11-08): the complete-near native800 branch
has reached t=1.0205 after a verified local seed-box search correction. The
same-time native400/800 comparison reduces common-fluid wake density
discrepancy against finite SU2 from1.258% to.701%, and actual excluded
wall-area excess from.02731745 to.01083429. No cell intrusion is found.
This does not establish whole-method convergence: the timestep caps also
differ (.00025/.0005), earlier temporal qualification was local, and a
1.683% finite-step volume-increment residual remains. Instantaneous RHS
consistency and energy-ledger closure are separate tests. Native800 is being
continued from its accepted1.0205 state to5; inspect spatial800_v1/native800_t5_v1
and its control/completion files for current status, not this launch note.
The complete-near native400t5 result, rather than older chart-only histories
below, supplies the current matched wake baseline. The original historical
figures and tests remain preserved.

Saved native400 and SU2 t5 velocities both have opposite-sign upper/lower
wake shear circulation. This supports the physical shear-layer interpretation
but is not a KH growth-rate measurement or proof that each density clump is
a vortex. Neither visual smoothness nor single-time density L1 determines
which solution better resolves KH. Evaluate temporal roll-up and resolution.

Current endpoint status (2026-11-08): qualified native400 completed t=5,
exit0, checkpoint `cylmat.020073.dat` in
`cylinder_late_failure_v1/chart_translation_v1/native400_t5`. Retained masses
are exact and U/V positive. Actual native wall-area excess at t5 decreases
3.43039% to1.19555% from200 to400; maximum radius .588751 to .529247.
SU2 h.025/dt.00025 now also completedt5 and passed saved-state, inner-solve,
restart/hash and zero-recovery checks. Against that reference, common-fluid
density relativeL1 decreases5.906% to3.329%, and wake density13.863% to10.526%.
Native400 pressure drag coefficient1.31882 versus SU2fine1.35474 (-2.65%).
Central threshold shockx1.688 versus reference1.690. The two SU2 spatial meshes
themselves differ2.38% in interior density and4.24% in wake density; these are
not rigorous error bars. Local nose Cp, approximate wall motion and wake
texture remain unresolved. Full comparisons and source provenance are in
`cylinder_endpoint_compare_v1`.

Integrating the fine reference over the SAME native cells does not remove
the wake discrepancy: native400 pointwise density difference10.43% becomes
9.50% for cell averages (selected full wake cells; region differs from above).
Thus polygon display/coarse averaging alone is not a sufficient explanation.
No smoothing or physical-state modification is applied by this diagnostic.

Whole native400 restart lineage, without discarded overlap, has stationary
wall physical relative energy defect -5.44838e-5 at t5, with ledger closure
-1.54e-14. Accepted-step timers total6919.918s, not full process time.
Matched-accuracy <=2xGIZMO remains unproven. Keep cylinder first; no3D yet.
Larger timestepcap.001 versus.00025 now completes saved.01-to5, including a
verified restart after an unexplained session143 termination. Full physical
relative energy defect is -5.91122e-5. It is not adopted: common-space wake
density changes6.49% between caps while both remain about10.5% from SU2fine.
Accepted-step cost is2102.800s, with a38.938s discarded tail recorded separately;
original process wall time is unavailable. This is not a GIZMO comparison.
Existing previous-CFL proposal reuse independently reduces a matched short
window's process time35.96→25.30s with rho/P differences around1e-9 and all
stage checks retained. Subsequent cylinder runners enable that existing feature,
not the larger timestep cap. Its whole-run/3D benefit remains unproven.

Frozen heat/work decomposition and matched t=.02 reference integration now
separate actual startup compression from unproven spurious-heating claims.
Selected native cells agree with reference means to.154% density/.617% pressure
at .02, despite substantial instantaneous wall-associated entropy generation.
The later low-density cohort mismatch must be localized at intermediate times.
See cylinder_endpoint_compare_v1/README.md for scope and raw-geometry limits.

### Historical progress records

Earlier snapshot (superseded by the Status section above):
original native400 stopped at4.05125. A centroid-roundoff-scale correction
passes the saved4.02125→4.06 window and leaves normal200/400 windows bitwise
unchanged. A later continuation stopped at4.0835 because periodic chart
translation loses coordinate low bits. A common C/CUDA translation correction
is implemented and being qualified; do not report native400t5 as complete.
SU2 h.05 both dt.0005/.00025 runs completedt5; density/pressure relative
time-step differences are9.80e-6/1.27e-5. The h.025 reference is still live.
Actual GPU wall segments confirm the coarse200 wall envelope exceeds the
specified circle's area by3.43% at t5, concentrated behind the cylinder.
This is a geometric boundary discrepancy, not just a density rendering mask.
Matched spatial convergence and morphological acceptance remain unproven.
Follow-up: chart correction passes36 CPU/CUDA vertex and24 CUDA clipping
cases, both normal native windows bitwise, the4.06→4.1 failure window, and
four existing2D short regressions. A new native400 continuation from4.1 is
live in cylinder_late_failure_v1/chart_translation_v1/native400_t5.
At matchedt1.06 the actual wall-area excess decreases7.14%→3.57% from200 to400;
stored and reconstructed volumes are bitwise equal. Both h/dt changed, and
the requiredt5 fine morphology comparison is still pending.

Not complete. Mach.1 mixed-boundary diagnostic reaches t=.01 on3912 cells;
its radial image boundary moves and has nonzero wall work. Native RunCylinder
material stage loop, long-time open-boundary events, Mach3 target and
independent morphological comparison remain unqualified.

After cylinder, KH/RT require GIZMO-level resolution and nonlinear evolution.
Sedov source currently combines a Cartesian particle array and four central
heated cells: separate particle-layout anisotropy from square deposition
before attributing all angular structure to the hydro solver. Only selected
validated comparison figures go in the manuscript; retain diagnostics outside.


## Updated evidence and reference scope

Historical baseline recovered and audited in cylinder_reference_audit_v1.
Published Guermond et al2018§5.7 has different domain/transverse walls and
Ut/D30, not the manuscript Ut/D5. Merely reaching manuscriptt5 will therefore
not prove published-level completion; use a matched independent comparison.
Finite-resolution radial-image wall errors may be accepted only with physical
error/refinement evidence, not by requiring every approximate power face to
be mathematically stationary. Physical wall contact impulses are distinct from
centroid steering, but need explicit event, impulse and energy accounting.

New actual Mach3 open-stream attempt: target5, accepted.081, stoppedbefore
committing wall penetration. Fourfold/sixteenfold dt reduction at the same
checkpoint still penetrates around.0813. Current status remains NOT COMPLETE.
Detailed data in cylinder_mach3_native_v1. Priority remains cylinder first.

Latest correction: the separate comparison runner had omitted the original
Exam/Cylinder/cylinder.c::cyl_postStage wall treatment. It now calls that actual
function at both SSPRK2 stages, with its kinetic-energy loss explicitly
ledgered. The old .081 penetration stop is passed. After an outlet-state fix
and a roundoff-certificate correction, the 200x40 comparison is continued from
the verified t=1.252 checkpoint toward5. This is not yet full RunCylinder,
the400x80 refinement, or a matched published morphology comparison. See the
last appendices of MATERIAL_RT_PORT_AUDIT_2026-11-06.md for exact provenance.


Latest endpoint: the material-departure boundary comparison reached t=5 on
200x40 sites. Outlet density stays at its IC value to5.35e-11; retained masses
are exact, U/V stay positive, and final generators are outside the obstacle.
Whole-run stationary-wall energy error is1.50e-4 relative, including image-wall
work and the original contact correction loss. Maximum transient relative
volume-increment residual .00861 still requires refinement. The current
periodic transverse domain is not an isolated-cylinder reference. A native
SoA helper now reproduces the same mixed RHS and short C/CUDA evolution;
actual RunCylinder stage/padding/stream/checkpoint wiring remains incomplete.
These are concrete advances, not completion of400x80/published morphology,
full native integration, multi-rank or <=2xGIZMO requirements.

Native wiring update (2026-11-08): RunCylinder now actually evolves the mixed
material boundary and restarts through owned-cell inflow/outflow. The t1-to1.06
comparison/restart differences are below2e-15 in densityL1; this is a port
qualification only. The original400x80 IC has reached.01 in a40-step native
pilot with fixed masses and positive U/V. No400x80t5 result or independent
morphology acceptance is claimed. The next cost task is the same anchored
weight solve on the existing CUDA PCG engine, keeping its physical operator
and true-residual criterion unchanged; then continue the saved refinement.
See cylinder_full_native_v1/README.md for the complete evidence and limits.

CUDA qualification update (2026-11-08): the same anchored FP64 weight solve is
now selected in the material profile after native200/400 time-evolution
comparisons and all four existing KH/RT/Sedov/Noh regressions. Density
differences are at roundoff. At400, median step cost fell from.512 to.286s on
the shared Ada, without changing physics or residual thresholds. A single
400x80 native continuation from t=.02 to5 is launched in
cylinder_cuda_anchored_v1/native_400_t5; inspect its live process/log or final
audit before reporting state. This does not yet satisfy endpoint, morphology,
wall/volume refinement, independent reference, or<=2xGIZMO requirements.


Transient geometry update (2026-11-08): native400 and the fine SU2 reference
have both reached t5; the historical launch status above is superseded.
However, final-only nonpenetration is insufficient. At native t=.27025 the
original wall admits .00288979927 area inside the cylinder (minimum radius
.47591206); centres remain outside. Completion now explicitly requires whole
cell nonpenetration at every geometry stage and saved intermediate times.
The identified cause is missing self images beyond the fixed4h padding range.
The C/CUDA support-completeness candidate passes a frozen repair and0-to.01
native regression, while .01-to.27025 time evolution is being qualified.
Neither stationary-wall accuracy nor whole-run morphology is yet accepted.
Full provenance is in cylinder_endpoint_compare_v1/README.md and wall_support_v1.

The support-completeness candidate has now completed t=.27025. Every evaluated
native geometry passed its tangent-halfspace certificate; an independent
whole-cell audit at the endpoint finds no disk intrusion. Stored and rebuilt
volumes match exactly. This closes the demonstrated early missing-image case,
not the full t5 morphology, stationary-wall work, periodic support or speed
requirements. Subsequent checkpoints must retain these geometric checks.


Temporal wall support update (2026-11-08): the self-only supplementation
candidate reached .5205 without circle intrusion but failed temporal geometry
consistency: active neighbour images were dropped, producing a dt-independent
.24% volume jump. The replacement includes the complete R<r<2R radial family
and certifies whole-cell support plus strict inactivity of the limiting
centre image at every stage. It removes adaptive supplementation/rebuilds.
The compatible .27025-to.28425 native refinement now reduces maximum relative
volume defects1.50e-5→4.16e-6→7.36e-7 asdt halves, with fixed masses, positive
states and independent nonpenetration. This closes the local deletion defect,
not the complete cylinder result. A .5205 continuation is being qualified.
See wall_support_v1/transition_v1/full_near_v1/production_v1 for actual status.

The local complete-near continuation has completed .5205 with no circle
intrusion and conserved retained masses. Matched-reference selected-cell
rho/P differences remain3.16%/4.41%, not improved overall accuracy. A later
contact-associated finite-step volume defect remains. Because the inherited
.27025 state passed through earlier candidate deletions, final method/history
qualification now starts from the compatible saved.01 state. The .27-to.52
local branch remains diagnostic. Four non-cylinder short regressions are
bitwise equal to the latest baseline; no3D or GIZMO completion is implied.

Earlier-history re-evolution has now established exact restart equivalence
through .27025 (all physical/cached values and exterior state equal; only
serialized linked-list addresses/checksums differ). The completed .5205
candidate is retained. Its single1.0205 continuation is live in production_v1;
check actual process/audit before reporting completion. Early contact-related
local volume defects remain; full cylinder and matched GIZMO gates are open.

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
