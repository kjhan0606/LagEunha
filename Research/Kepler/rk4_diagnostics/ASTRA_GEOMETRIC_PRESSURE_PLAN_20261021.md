# Astra corrected plan: geometric volume / compatible pressure and heat

2026-10-21. ONE end-of-shared-factory bundle evaluation with a mathematical
consistency correction. Driver codes/runs. No new microreview is needed
during implementation of this coherent bundle.

## Evidence and withdrawn first formula

Shared-state IC/CPU/CUDA/checkpoint and sparse/global work gates pass;
max chi=.003078319 fails unchanged.002. Saved load decomposition attributes
.99999999165 of worst divergence error to the nearly constant pressure
weak-load residue. Gravity remainder chi5.50001e-5 passes. No significant
raw-Gram normalization discrepancy was found.

Astra initially proposed fixed Pref=P0 residual force/work plus boundary
geometry. Driver identified that this only guarantees nullspace for
P=Pref, not EVERY constant P. For P*=c !=Pref, the nonzero defect
(c-Pref)*(Q[cofF gradN]-G) remains. Astra WITHDREW that claim and the old
fixed-Pref thermal potential as the implementation contract. Do not
implement force-only subtraction or omit a pressure-projector time term
while claiming the old potential unchanged.

## Final selected spatial discretization

For each physical cell, SAME positive reference quadrature weights w_q:

    A=sum w_q, alpha_q=w_q/A (fixed projector weights),
    Pi(P)=sum alpha_q P_q,
    DeltaV=V_boundary - sum w_q J_q,
    nu_q=w_q J_q + alpha_q DeltaV.

Then sum nu_q=V_boundary. Require EVERY nu_q>0; never clip or redistribute
on failure. This explicitly changes the THERMODYNAMIC DISCRETE VOLUME.
It is not merely a different output ledger. Validate nu versus actual
w_q J_q with existing physical accuracy limits.

Keep SAME saved x/v, material masses m_q, PHYSICAL U_q, gravity, reference
mesh, quadrature, sparse mass matrix and CG. The previous control stores
thermal state implicitly as the original adiabatic control law; export its
unchanged old point U_q ONCE before the new EOS, rather than reinitializing
U from nu or using the control law as a production pressure reset.

Define physical thermodynamics on the selected discrete volume:

    rho_q=m_q/nu_q,
    P_q=(gamma-1)*U_q/nu_q,
    F_a^P=sum P_q * partial(nu_q)/partial(x_a),
    Udot_q=-P_q * nudot_q.

Equivalent pressure projector form:

    F_a^P=Q[(P-Pi(P))*cofF gradN_a]+Pi(P)*G_a,
    sum Udot=-Q[(P-Pi(P))*J divv]-Pi(P)*Vdot_boundary.

Pref may remain FIXED at original P0 for stable arithmetic splitting:
Pi(P)=Pref+Pi(P-Pref). It is not adjusted with error/configuration.
Nullspace follows from the fixed projector/nu, not Pref.

For EVERY spatial constant P=c, F=cG and Udot=-cVdot. External pressure
and internal paired faces must use the SAME geometric derivatives, so the
constant-pressure cancellation is algebraic. Do not threshold P-Pi(P).

At fixed entropy, U_q=A_q nu_q^(1-gamma) differentiates to the stated
force/work. There is no omitted dot(Pi P)*DeltaV term: do not invent
-(Pi P)*DeltaV as an extra independent thermal potential. This entropy
relation must not be used to overwrite/reset the saved current U_q.
Shock entropy and long evolution remain separate later obligations.

## Boundary volume and its ACTUAL derivative

Define V_i^b from canonical oriented physical cell boundary edges, evaluated
once and scattered with paired signs. Same-cell internal triangulation
edges cancel topologically. A nonpolynomial edge can be written

    x(t)=a+t(b-a)+delta(t), delta(0)=delta(1)=0,
    V_e=.5 cross(a,b)
        +Q_e[cross(delta,b-a)+.5 cross(delta,delta')].

Use translation/rotation-invariant cell expressions and stable differences.
G_ia=partial(V_i^b)/partial(x_a), and Vdot_i^b=sum G_ia dot v_a.
DIFFERENTIATE THE CHOSEN DISCRETE FUNCTIONAL, including both delta and
delta' terms. A separately quadrature-integrated N_a cofF n is not
automatically its derivative. Compare that cofactor integral independently
and account for the difference as geometry error.

Internal face baseline force/work/torque is one canonical record with
opposite scatter. External pressure baseline must use the same G; do not
add a separately approximated -P0 n that breaks cancellation again.
Variable external residual force and its work/torque must share a native
operator, with physical-traction error checked independently.

## Same-state implementation/verification bundle

1. Read SAVED initial and finite shared controls. Do not construct a new
   trajectory, replace old acceleration, refit states or change degree /
   mesh / quadrature / precision / criteria.
2. Freeze old point masses/U; compute Vb/G/Vbdot, DeltaV, nu/EOS pressure
   and fixed pressure projector. Save all physical/geometric differences.
3. New pressure RHS and SAME sparse mass solve. Preserve old a as a control.
4. CPU/OpenMP/MPI/CUDA checks at both saved states; arbitrary constant
   pressure, gravity-off/linear/Plummer conditions and original chi gates.
5. Complete cell-wise numerical traction, energy exchange, actual internal
   angular-momentum/rotation rates and positive physical U/J/nu. No
   independent spin/kinetic reservoir, viscosity, heat top-up or steering.
6. Keep checkpoint epoch/owner/caller preservation and original state SHA.

Additional mandatory gates: sum nu=Vb, every nu positive, nu-wJ error,
new EOS versus original MATERIAL pressure error, Vb versus actual volume,
unchanged physical U integrals and old curvature (max.002/every mean.001).
Compatible work alone is insufficient. On any failed accuracy/positivity/
geometry contract, preserve witness; do not apply only force or correct U.

On pass, next is SAME target68526 complete transfer with prospective shared
state. The earlier inverse-domain three failures remain unclosed until
actually verified. Original CUDA/MPI long t6/8/10 and GIZMO longest orbit
goal remains ACTIVE and UNACHIEVED.
