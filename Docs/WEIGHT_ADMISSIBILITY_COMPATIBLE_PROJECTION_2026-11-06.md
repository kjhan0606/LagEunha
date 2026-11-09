# Candidate admissibility projection: experimental native pilot, not adopted

The seeded material KH stop shows that prescribed volume matching does not preserve the native admissible weight set. Existing margin beta=.9 is retained. No weight clipping, generator steering, changed viscosity or new fit coefficient is authorized by this note.

At a fixed state write the current compatible semidiscrete operator as

    Vdot = C v + D P
    F = C^T P - K v
    H_total = P^T D P + v^T K v
    Udot_total = -P^T Vdot + H_total.

Here D=Jp^T diag(A/(ZL+ZR)) Jp and K=Jv^T diag(A Zeff) Jv are symmetric positive semidefinite, using the same reconstruction and transpose. C represents the acoustic weighted mean-velocity owner-volume map. These are algebraic operators, not a proposal to assemble dense production matrices. L=dV/dW and B=dV/dx are the geometric maps, so L Wdot=Vdot-Bv.

For a single active pair inequality s(Wi-Wj)-beta|xi-xj|^2<=0, define a=s(ei-ej), and b^T v=2 beta (xi-xj).(vi-vj). Its tangency condition is a^T Wdot<=b^T v. On connected geometry solve Lg=a with zero mean, giving

    g^T Vdot <= h^T v,  h=B^T g+b.

Set c=C^T g-h, S=g^T Dg. If S>0 and z=c^T v+g^T DP>0 at the boundary, lambda=z/S. A compatible projected operator is

    C' = C-(Dg)c^T/S
    D' = D-(Dg)(Dg)^T/S
    Vdot' = C'v+D'P
    F' = C'^T P-Kv
    H' = P^T D'P+v^T Kv.

Then g^T Vdot'=h^T v, total energy balances exactly at the semidiscrete level, and D' is positive semidefinite by the D-inner-product Cauchy-Schwarz inequality. If inactive, use the original operator. Translation invariance is retained when the original C annihilates uniform translation and the geometric boundary functional has the same property; global volume cancellation follows ones^T C=ones^T D=0. These conditions must be checked on actual native operators.

Interpretation: this changes the compatible pressure response and pressure-jump relaxation when geometric freedom is exhausted. It is not enough to alter Wdot alone. It is a constrained numerical closure, NOT established physical truth. Whether its constraint force respects the user's intended no-steering physical interpretation and yields accurate hydro remains unproven; no production force changed here.

## Completed check and limits

Research/MaterialVolume_20261105/weight_constraint_projection_20261106_v1/audit.py tests260 active random connected graph cases at n=4,8,16. Tangency, momentum, volume and energy identities hold to floating-point error; D' smallest eigenvalue>=-3.8e-15, all tested heat>=0. Naive volume-only modification gives large energy defects. Python is algebra verification, not simulation.

Still required before any native trajectory use:
- Audit actual saved-state C,D,B,L and force/heat identities, including nullspaces.
- Handle S=0 or tiny without an arbitrary floor; infeasible constraints must fail explicitly.
- Multiple simultaneous constraints with complementarity; single-pair algebra does not solve the observed multi-hit case.
- Define activation/event handling and finite-step geometric admissibility under SSPRK2, including changing neighbor sets. Continuous tangency alone is insufficient.
- Construct a justified local heat allocation preserving the total identity and positivity; the global quadratic proof alone does not supply positive per-cell heating.
- Bound runtime: extra L solves and active-set work may defeat the <=2x GIZMO goal. No performance or dimensional qualification yet.

This is a concrete candidate to evaluate against the saved stop state, not a reason to alter the current margin or declare KH solved. All native trajectories and checkpoints remain unchanged.

## Native pilot update

The earlier algebra-only status above is historical. Native C projection with CUDA geometry has now been tested on the seeded KH saved stop, including simultaneous constraints and positive local heat. The latest constraint-union continuation reaches t=1.28 at two time steps. See `Research/MaterialVolume_20261105/native_constraint_union_20261106_v1/README.md` and `comparison.json`. Local pressure changes by about 10.3% over t=1.271825–1.28 while halving dt changes pressure by at most 0.132%. Cost is about .25 seconds per step. Neither physical accuracy nor the runtime target is established; no common-method adoption, full GPU/MPI qualification, RT or 3D completion is claimed. General complementarity and finite-step admissibility remain incomplete.
