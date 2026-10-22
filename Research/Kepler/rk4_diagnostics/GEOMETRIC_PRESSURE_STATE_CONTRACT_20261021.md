# Compatible geometric pressure: saved-state contract

Driver implementation of the corrected Astra bundle, 2026-10-21.
Authoritative code identity remains LagEunha HEAD
8e25dfe0c2e41a7374c11e2e3102da2a8b028625. These are separate native
experimental controls; original repository code and retained states persist.

## Thermodynamics and actual force

For each physical cell, fixed positive reference weights define
alpha=w/sum(w). Canonical boundary volume Vb and its actual discrete
configuration derivative G define

    nu_q=w_q*(J_q+(Vb-QJ)/sum(w)),
    P_q=(gamma-1)*U_q/nu_q,
    Fp=Q[(P-Pi(P))*cof(F)*grad(N)]+Pi(P)*G,
    Udot_q=-P_q*w_q*(J_q*div(v)+(Vbdot-Q[J div(v)])/sum(w)).

Pi(P)=sum(alpha P). The reference pressure P0 is fixed and used only for
stable subtraction. Every nu must remain positive. The previous physical
U, reference density/weight and x/v are exported before this EOS change.
U is never recomputed from nu or reset after a force update. No trajectory
is advanced in this bundle. Separate manufactured constant-pressure unit
tests change their own input U only; they do not change either real state.

The same positive 64-angular/28-radial volume rule, 64-point canonical edge
rule, reference basis, sparse mass, original PCG accuracy and native gfs
stress-work operator are retained. CUDA actually evaluates pressure/EOS,
Piola stress, gravity and thermal work, with the same shared coefficient
and reference tables as CPU/OpenMP/MPI. Canonical boundary geometry is
computed once per reference edge on the host; this is not a claim that
all geometry or the sparse solve runs on CUDA.

## Open-edge derivative and face traction

For x(t)=a+t(b-a)+delta(t), delta endpoints zero, chosen finite-rule volume
is .5 cross(a,b)+Q[cross(delta,b-a)+.5 cross(delta,delta')]. Both delta and
delta' terms are differentiated. The derivative of an OPEN edge's area
contains endpoint gauge terms; it cannot be used directly as physical
face traction. Use

    H_A=Gedge_A-.5*Rcw(a), H_B=Gedge_B+.5*Rcw(b),
    H_other=Gedge_other.

The endpoint terms cancel around each closed cell, so the signed sum of
canonical H is exactly the cell volume derivative G. Each paired H is
scattered with opposite signs. H is also compared independently with
Q[N Rcw(x')]; no assumed quadrature integration-by-parts identity.
Exterior pressure retains the original material face pressure P0*J^-gamma.
Its fixed edge-weight mean multiplies H, and its residual multiplies the
canonical point traction. Work and torque use these actual same loads.
Translation acts on P3 nodal coefficients only, not enrichment coefficients.
The preliminary v1 translation diagnostic erroneously summed both; v2
corrects that diagnostic and adds it to the gate. No physical kernel change
was needed for that correction.

## Resolved rotation and hybrid exchange

Cell reaction is Mcell*a minus the actual cell pressure/gravity/exterior
load. Global shared reactions must cancel; private zero-trace bubble
reactions must vanish. Reaction work and torque are paired with the SAME
stored shared v and x coefficients. No fitted force or heat is applied.

Actual mass-centred spin S=L-Rcm cross P, inertia I=int|x-Rcm|^2 dm and
their rates are computed from the represented x/v/a. omega=S/I and
Erot=S^2/(2I) are derived portions of the existing kinetic energy, with
Erot_dot=omega*Sdot-.5*omega^2*Idot. They are not an independent spin state
or extra energy reservoir. The remaining internal kinetic-energy rate is
reported explicitly. Kepler Re=infinity; no viscosity, damping or steering.

## Native gates and limits

Unchanged initial curvature: maximum chi<=.002, every cell mean<=.001.
Energy-rate residual/U<1e-8 globally and per cell; angular residual/2K<1e-11;
CPU/CUDA force difference<1e-11. Positive J/nu/U and untouched x/v/U are
mandatory. New physical volume/EOS differences must be<1e-10; boundary
functional/cofactor and external traction differences<1e-8. These additional
accuracy bounds do not replace or relax the older curvature/work gates.
Canonical translation and rotation identity defects must be<1e-12.
Gravity-off and linear gravity are same-state operator controls only.

Restart tests still preserve caller on wrong epoch/duplicate/missing ID.
The frozen per-point physical table is a separate, epoch/versioned thermal
companion (triangle,X,Y,weight,rho,U), protected with input/output SHA.
The current shared-state checkpoint by itself is not a future complete
production thermal restart. Physical evolution needs persistent point U
and its time integrator; it must not reuse the old control's constitutive
formula to regenerate thermal state.

These are native engineering accuracy checks, not outward interval bounds
for analytic caches/functions or a proof of whole-target transfer. Exact
C0 references and compiler obligations in SHARED_TRACE_STATE_CONTRACT
still apply. Full target68526, inverse-domain failures118/246/254,
production timestep/restart, historical t6/8/10 and GIZMO longest-duration
validation remain subsequent work. No long-Kepler stability claim follows
from these instantaneous finite controls.
