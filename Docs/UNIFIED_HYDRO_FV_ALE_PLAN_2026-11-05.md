# Common conservative FV--ALE candidate: grouped Astra plan

The full user objective is unchanged: explain the 1D contact error, use one
method in 1D/2D/3D, retain self-consistency and simplicity, rerun the 2D suite,
and measure wall time at most twice matched GIZMO. This plan follows the
specific unlimited kinetic-recovery rejection in `90679f2`. It does not
promote a new solver or replace any accepted trajectory or manuscript result.
One grouped Astra planning call selected the following candidate and tests.
The driver implements it; no per-test planning calls are required.

## State and physical velocity

Use `d+2` extensive Euler variables `Q=(m,p[0..d-1],E)` on each actual power
cell. Define `rho=m/V`, `ubar=p/m`, and

    P=(gamma-1)[E-|p|^2/(2m)]/V.

This is an explicitly declared finite-volume mean-state approximation, not
an exact recovery of integrated physical thermal energy from arbitrary
subcell data. The measured h^2 variance-transport term remains an error to
quantify in actual cold flows. No kinetic deconvolution, extra thermal state,
weight energy, or reservoir is added.

For a NEW analytic IC, use actual-cell mass/momentum integrals and set
`E=|p|^2/(2m)+integral P/(gamma-1)`. Record the difference from the continuum
total-energy integral as initial spatial projection error. Do not convert
old accepted checkpoints or reset E after evolution.

Use centroid displacements for the common LSQ reconstruction. Move generators
with the reconstructed physical velocity `qdot=ubar+G(q-c)`, with the same
fixed reconstruction limiter as the face velocity. This adds no steering.
The reconstruction/limiter still requires an explicit implementation and its
own dimensional checks; it is not established by the fan audit.

## Complete actual-face flux

Freeze the reconstructed lab-frame HLLC/HLL fan during each weight-rate solve.
Sample `F_ALE(w)=F_k-w U_k` for mass, every momentum component, and energy.
Never discard face mass transport or the numerical HLL residual. Keep the
original Davis HLLC primary bounds and star-state admissibility check. For
fallback use the COMPLETE HLL fan with normal-mean-frame outer bounds:

    b=(uLn+uRn)/2.
    a=max(|uLn-b|+cL, |uRn-b|+cR).
    sL=b-a, sR=b+a.

Limit rho/P slopes for positive states at the actual quadrature points;
cell Q is unchanged. Fan admissibility is not a theorem for a full timestep.
Fix face integration to a point in 1D, three-point Gauss in 2D, and positive
degree-four-exact Duffy triangle quadrature in 3D. The existing 3x3 Duffy map
may be reused after verifying the required degree. The polynomial-flux
motivation is [Zier--Springel section 4.2](https://arxiv.org/html/2205.07916#S4.SS2),
not a certificate of thermal accuracy or code performance.

## Weight-linked mass without a hard constraint

At birth, store positive, immutable `mu_i=m_i0/V_i0^(2/d)`, `m_i0`, and `W_i0`
in restartable state. Set `z=Wdot` by

    R_i(z)=mu_i z_i-sum_faces integral Fmass(wq+(z_i-z_j)/(2r))=0.

The Jacobian is `diag(mu)+Lrho`. The native graph solve must stay in one
backend, rather than performing GPU--host JVP copies each PCG iteration.
Freeze caps at eight Newton iterations and 64 PCG iterations per Newton
step. Do not change mu, increase caps, or subtract an RHS mean to force a
pass. Do not recenter W during evolution.

Mass is updated by the complete face flux. The law gives `mdot=-mu*z`, so
`m+mu*(W-W0)=m0` is an invariant to the integrated true residual. Log its
actual defect, without overwriting mass from the invariant. Maintained
`z=0` implies fixed mass; maintained `W=0` is Voronoi. Forcing z to zero in
arbitrary flow while retaining unconstrained ALE mass flux is a different
model and does not satisfy this reduction.

## Time integration and geometry

Use SSPRK2 for Q, q, W, reconstructing actual geometry at each required
stage and final state. Reuse final geometry on the next step. Fix initial
wave CFL to `0.4 min_i 2V_i/sum_faces(A*a_face)` with
`a_face=max(|sL-wn|,|sR-wn|)`. On every forward-Euler ray, bound time by the
FIRST zero of `m(s)>0` and `2m(s)E(s)-|p(s)|^2>0`; use a fixed safety factor
(implementation registration must specify it before PDE execution).
Validate the actual trial state and V, including the second stage and final
combination. Reject the entire step and halve dt upon invalid state or
geometry, with a cap of eight retries. No floors or energy repairs.

Record

    deltaV=V(q1,W1)-V(q0,W0)-dt*(Vdot0+Vdot1)/2.

SSPRK2 on q,Q,W does not automatically satisfy exact geometric conservation.
Measure its time order with one registered dt-halving control; do not claim
exact free-stream preservation. For cold flow also disclose cumulative
`sum |E deltaV/V|` relative to thermal energy. This is a diagnostic budget,
not a rigorous pressure-error bound.

## Ordered native evolution and fixed stopping criteria

The same settings, shared dimensional functions, and native binary path must
be used for the corresponding tests. Initial dimensions differ only through
the actual geometry adapter and component count.

1. Fan/boost/conservation and stationary/moving contact, including small
   actual native contact in all three dimensions.
2. Existing nine 1D tests at N200 and their full original endpoints; Sod,
   strong collision and Shu at N400 as well. Shu uses the corrected original
   `sin(5*pi*x)` condition and actual-interval independent reference.
3. Actual 2D axis Sod 64^2/128^2 and 3D axis Sod 8^3/16^3, all through t=0.2.
4. Actual Gresho 64^2/128^2 through t=1, and one dt/2 control at 64^2.
5. Periodic sine shear 64^2/128^2 through t=0.25, rho=1, P=1e-6.
6. If these pass, KH 64x128 through t=2 with the same candidate.

Every accepted state requires positive m,V,P, closed mass/momentum/energy
including boundary ledgers, and the mass--weight invariant within measured
residual accumulation. Noise in translating regular contact, nondecreasing
shock-refinement errors, or solver/retry caps fail the candidate. Compare
Gresho using physical reconstructed velocity at the generator AND declared
cell-mean norms; do not compare old nodal velocities to new ubar. Gresho must
improve with resolution and not worsen the existing baseline. Cold sine
shear targets volume-weighted pressure relative L1 <=1%, positivity and
resolution improvement; time error dominating spatial error is unresolved,
not passed. No test-specific mu, CFL, quadrature or limiter adjustments.

Do not graft unqualified gravity onto this FV branch or restart old 20pi
states under a different state interpretation. Long Kepler, the remaining
2D suite, MPI and 3D GPU qualification remain later requirements.

## Measured performance, not projected GPU speed

After physics gates, sequentially compare Gresho and cold shear with GIZMO
at matched analytic IC, EOS, domain, output, resolution, endpoint and
resources. Report initialization plus evolution and evolution-only times.
Split native time into geometry, reconstruction/quadrature, fan, weight
solve and rejected work. Any measured wall time above 2x GIZMO misses the
user target. A CPU prototype is not qualified by speculative GPU speedup.

## Current implementation boundary

`Exam/gfs_fv_ale.h` begins the dimension-common primitive/frozen fan and
forward-ray cone calculation. Its fixed five-slot payload embeds d+2
independent variables with zero unused momenta; energy is always slot four.
`fv_flux_audit.c` checks these native pieces. Optimized native C/C++ and
actual Ada audits pass in `UNIFIED_HYDRO_FV_PRIMITIVES_2026-11-05.md`.
No production geometry, weight solve, time update, restart ABI, or FV
trajectory is enabled by that header. The separate Shu IC repair retains
the existing W0 variational operator and implicit KDK; it is not a result
for the new FV candidate. The complete objective is unfinished.
