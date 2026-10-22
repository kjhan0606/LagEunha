# Compatible pressure/heat bundle: native saved-state results

2026-10-21. Approved corrected Astra proposal implemented by the driver.
Read GEOMETRIC_PRESSURE_STATE_CONTRACT_20261021.md and machine ledger
GEOMETRIC_PRESSURE_COMPLETE_LEDGER_20261021.json for definitions and scope.

## What changed

The physical point mass and U are retained. A fixed positive reference
projector defines nu=wJ+alpha*(Vb-QJ); EOS, pressure force and thermal work
are derived together from this same discrete volume. Canonical boundary
volume is differentiated, with endpoint gauges removed for physical face
loads and opposite signs on paired faces. The original native gfs stress
operator is reused in actual CPU/OpenMP/CUDA kernels. No damping, heat
repair, viscosity, independent spin reservoir, steering or trajectory reset.

Resolved mass-centred internal angular momentum, inertia and their rates,
derived rotational kinetic energy/rate and remaining internal kinetic rate
are now recorded per cell. Hybrid reactions obey the same stored x/v work
and torque. These are cell numerical exchange ledgers; they are not a
reconstruction of independent physical traction functions on every face.

## Fixed native controls and results

Two SAME SAVED states from shared_trace_rhs_v1_mpi4: initial_state.bin and
shared_state.bin.24 physical cells/384 source triangles/817152 original
positive volume points/3797 scalar shared DOFs. These are finite controls,
not a full128^2/256^2 evolution. Same volume/edge rules, mass matrix and
PCG accuracy, no mesh/order/precision/tolerance scan.

geometric_pressure_v2_mpi4 and geometric_pressure_v2_mpi1 both exit0.
MPI point ownership differs while original frozen owner4 metadata is
preserved; these are read-only state/RHS comparisons, not a general
MPI repartitioned evolution/restart demonstration.

| Check | Result (worst initial/finite unless specified) |
|---|---:|
| Previous initial max chi, original nearly constant-pressure load |.003078319055743 FAIL|
| New initial material chi |5.50000996964e-5 PASS(.002)|
| New initial every-cell mean chi max |1.50359602939e-6 PASS(.001)|
| Physical point U change |0, exact native equality|
| Thermodynamic volume relative change |4.43764e-16|
| EOS pressure relative change |8.47033e-16|
| Global energy-rate residual /U |1.69193e-15|
| Every-cell energy-rate residual /U |6.54530e-14|
| Global angular-rate residual /2K |1.48577e-17|
| Every-cell angular-rate residual /2K |1.24995e-16|
| CPU/CUDA force difference |5.13228e-13|
| Actual mass-equation relative residual |9.29915e-17|
| Global hybrid-reaction residual /load scale |8.84268e-17|
| Private zero-trace bubble reaction /load scale |8.84268e-17|
| Canonical derivative/cofactor difference |8.13259e-17|
| Exterior traction versus original Piola control |2.08784e-16|
| Native arbitrary constant-pressure residual force |1.62758e-16 relative|

Every sampled J/nu/physical U positive. Small quadrature weights explain
small minimum point U; no floor or redistribution is applied.136 physical
boundary segments:56 paired,80 exterior;480 same-cell material edges cancel
topologically. Canonical traces/zero bubbles and checkpoint roundtrip,
wrong epoch/duplicate/missing caller-preservation pass.

Actual EOS/native stress/CUDA/global assembly were exercised for constants
.37P0,7.25P0 and1, including values different from the fixed P0 reference.
The algebraic canonical-face jump is zero; the real native full residual
above includes arithmetic/representation roundoff. Gravity-off, linear
and Plummer ledgers pass at both real saved states.

MPI1/4 frozen reference rows (mass inputs and physical U) are bitwise
identical per source. Saved x/v/basis/owner checkpoints equal the original
input bytes. New accelerations and loads have exactly equal CSV numerical
values, including all3797 global IDs. Small differences in long-double
postprocessing reductions are recorded separately. The original initial
acceleration was read from its SAVED CSV, checked against the unchanged
old mass/load equation, and retained in old_new_loads.csv. Finite original
pressure acceleration is a separate SAME-state linear solve, no new state.

Native bundle wall times5.30s(MPI4×OMP4)/13.45s(MPI1×OMP4), compile excluded,
on authorized lageunha/Ada. Those numbers are not full-grid evolution speed.

## Astra end evaluation: actual nu curvature was also required

Astra accepted the nullspace/H construction but identified a missing
obligation: the new EOS uses nu, so inspect its actual second time
derivative. Implemented separately, without rerunning the mass/RHS solve:

    Pddot/P=gamma*(gamma+1)*(nudot/nu)^2-gamma*nuddot/nu,
    nuddot=w*(Jddot+(Vbddot-QJddot)/sumw),
    Vbddot=G dot a+v^T Hessian(Vb) v.

The SAME chosen volume functional is quadratic; Hessian(v,v)=2*Vb(v) is
calculated directly at unchanged reference functions. No finite differences
or guard/precision scan. Saved initial/finite x/v and SAVED new acceleration
are read directly. geometric_curvature_v1_mpi{1,4} both exit0:

| Actual temporal quantity | Initial | Previous finite control |
|---|---:|---:|
| max absolute Pddot/P divided by omega^2 |5.50000996938e-5|.00128370035801|
| every-cell mean max of that quantity |1.50359602939e-6|.00107066875228|
| max nu-versus-material-J curvature correction |2.86249e-15|1.17388e-14|

Initial actual-nu curvature passes unchanged max.002/every mean.001.
The finite state is the OLD kinematic control x=x0+h v0+.5h^2 a_old,
v=v0+h a_old (h2e-4), not an accepted RK4 step or an exact circular Kepler
trajectory. Its nonzero absolute curvature is reported, including mean
.00107067 ABOVE the initial IC bound. We DO NOT claim that finite state
passes ideal-Kepler evolution accuracy. Only the additional nu-versus-J
curvature correction is tested there. The finite result.json zero chi in
the earlier pressure-only ledger was unmeasured, now replaced with null
in the machine audit; actual temporal values are provided above.

Worst initial actual-curvature point is source377/cell68600,
X1.01619449e-6/Y.0351448114. This is different from the removed pressure
residue witness at source5 in the exterior atmosphere. Neither witness
alone diagnoses the historical t6-10 inner-disk collapse.

A preliminary inline postprocessing assertion incorrectly applied the force
relative tolerance1e-11 to the near-zero nu-minus-J curvature correction.
It failed on relative2.86153e-11; the absolute MPI difference was only
3.35908e-25 in chi units. No native physical gate was relaxed or simulation
rerun. The full ledger preserves this diagnostic failure and both absolute
and relative differences. Existing chi/energy/positivity criteria all remain.

## Provenance / next connection

Source/build/binary/output/input SHA recorded; all original states and old
result directories preserved. v1 was preliminary. Its translation diagnostic
summed enrichment coefficients as if they translated; v2 restricts it to
P3 nodal translations and adds the identity gate. Physics kernel unchanged.
Frozen runners never edited while running; fresh directories refuse overwrite.
No authoritative repository modification, commit, or push in this bundle.
All native handles90245/42750/16281/33998/60880/20383 CLOSED; no jobs or
monitor loops left by this bundle.

Current result: native saved-state constant-pressure/heat compatibility,
initial actual-nu curvature and resolved cell work/rotation-rate gates pass.
This is NOT full long Kepler stability, analytic interval certification,
whole target transfer, a complete general thermal restart or production
integration. CPU boundary geometry/host sparse solve are not claimed CUDA.

Continue with ASTRA_SHARED_TARGET_CONNECTION_PLAN_20261021.md: native shared
canonical function + immutable thermal off-node function, whole target68526
all384 sources/all12 moments, no old broken-packet re-integration. Preserve
rho0/J (geometric push-forward) versus rho0/(J+delta) (EOS) distinction and
cache/reference error bounds. Old inverse-domain118/246/254 failures still
require real closure. Then full targets/cold reexpression/actual native
RK4 and historical t6/8/10/GIZMO longest-duration validation. Goal ACTIVE,
UNACHIEVED; no user approval or external blocker is pending.
