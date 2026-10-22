# Native moving-cell Euler RHS and screened weight-rate solve

2026-10-21. Previous goal turn was progress: enriched initialization resolved
the cold initial-pressure projection obstacle. This turn implements the actual
initial moving-cell Euler moment RHS and a coupled mass/weight RATE solve.
It also identifies a new, independently reproduced thermal tangent failure.
**Do not evolve this fixed-reference enriched candidate as a physics pass.**
No enriched timestep, dynamic-weight geometry or new restart has been made.
t6-10 and the longest GIZMO interval remain open.

## Identity and implementation

Native repository remains the worktree in `ENRICHED_EULER_RESULTS.md`, HEAD
8e25dfe plus uncommitted driver changes. Frozen `enriched_v5_build` binary SHA:
`cc0f82f8849d87c0e842f024169a90b02936d8a456213914fd03e693b5fe6f99`.
All native calculations ran on authorized lageunha CPU/Ada. Python launches
and analyzes; actual solver and independent reference calculations are C++/
CUDA. Existing sources, inputs, failure outputs and restarts are retained.
No commit/push or global-context modification occurred.

New native files:

- `Exam/gfs_ale_fan.h`: lab-frame HLLC fan, complete HLL fallback, all4 Euler
  components sampled with F_ALE=F_region-w_n U_region.
- `gfs_enriched_coupling{.h,_gpu.inc,_native.inc}`: canonical native support
  faces, same-point left/right Euler state and moment flux; full Euler volume
  flux/gravity; screened weight-rate Newton/CG; MPI reduction; separate
  geometric moment transport. Opt-in `GFS_ENRICH_COUPLED_RHS=1` inside the
  projection-only mode. No clipping, floor, steering, reference reset or
  equilibrium residual cancellation.

Generator velocity for this new operator is reconstructed physical fluid
velocity at the generator. Old point velocities are not silently reused.
MPI allgathers the initial basis/state once and uses CUDA local face work plus
global vector reductions; it is an audit implementation, not the final scalable
halo/migration runtime. Legacy default path passes its3-step replay bitwise,
including metadata (`enriched_legacy_v5/manifest.json`,5.104s).

## Stage mass/weight relation and Jacobian

Retain the relation selected in `ASTRA_AFFINE_MOMENT_PLAN.md` after replacing
its rejected thermal closure:

    m_i + mu_i W_i = m_i0,  mu_i = m_i0 / V_i0 > 0,
    z_i = Wdot_i,  R_i(z) = mdot_i(z) + mu_i z_i = 0.

For each canonical positive quadrature measure A_q, generator distance d,
normal n from i to j, face point x_f and physical generator velocities v_i,v_j,

    w_n = v_j.n - (v_j-v_i).(x_f-x_i)/d + (z_i-z_j)/(2d),
    mdot_i = -sum_q A_q Fmass(w_n).

All states and geometry are fixed inside this stage solve. A fixed Riemann fan
has dFmass/dw_n=-rho_region<0, also in the HLL middle region. Each point adds
a[[1,-1],[-1,1]] to the Jacobian, a=A_q rho_region/(2d)>0; complete J is
diag(mu)+L_a, symmetric positive definite. A global HLL fallback is selected
if either HLLC star state is inadmissible. All Euler components retain that
same fan; no pressure-only replacement. This is a mathematical stage statement
for these fixed reconstructions, not a dynamic tessellation existence proof.

Canonical fluxes conserve total mass, so the root also gives sum mu_i z_i=0.
W=0 gives m=m0 in the algebraic relation. The new dynamic path generally has
nonzero Wdot already at its initial W=0 state. It is not a proof that the
full Euler mass flux can be discarded on an arbitrarily fixed-W0 trajectory.
The legacy fixed-W0/fixed-mass path remains separate and unchanged.

`fan_contract.json`:15030 samples, modes1-5 including9 HLL fans; no monotonicity
violation, translated-flux discrepancy6.32e-13, flux continuity6.23e-13,
boost discrepancy6.08e-12. Bounded8-node graph Jacobian finite-difference
error5.61e-10, symmetry0, Rayleigh lower-bound ratio>=1.0414. Newton2 steps,
residual4.53e-16. No hard fixed-mass DAE, ridge or hidden mass/weight clamp.

## Actual native initial RHS results

Seven explicit runs:4096 uniform, rigid, boosted shear, varying-rho isobaric;
71080 Kepler MPI1, MPI4, and uniform velocity perturbation(.013,-.007).
All input hashes unchanged and all calls exited without time advancement.
Full candidate mode still uses original smooth mathematical IC, epsilon=.047
and the existing refined layout. No Hopkins/GIZMO profile was substituted.

Actual Kepler:142597 canonical neighbor pairs,2277312 positive face points.
Weight solve uses Newton1, CG34, no backtracking. Mean mass/momentum/gas-E
rates are canonical full Euler flux plus actual gravity; higher modes use
the same numerical flux at the same points. Spatial momentum modes carry
intrinsic angular information, with no scalar torque-power heat subtraction.

| Actual Kepler quantity | Result |
|---|---:|
| max abs(mdot+mu Wdot)/mu | 3.92e-16 |
| total mass rate | -5.66e-20 |
| weighted gauge rate | 2.02e-18 |
| sampled Jacobian symmetry discrepancy | 2.26e-16 |
| Rayleigh / diagonal-mu lower bound | 1.29531425 |
| independent secant / same bound | 1.29531425 |
| max abs(Wdot) | .00106389423 |
| MPI1/4 max Wdot difference | 2.06e-18 |
| MPI1/4 max normalized RHS difference | 4.58e-16 |
| velocity-boost mean-rate transformation discrepancy | 9.34e-15 |

Uniform moving-state GCL: normalized RHS minus geometric transport<=2.534e-13;
mean mass GCL2.10e-13. Actual rigid/shear periodic controls include2160/912 HLL
points and converge with the same screened solve. Their wrapped affine fields
have boundary discontinuities; global stationary-flow claims are not made.

Native wall times including startup, projection, global initial-state gather,
all RHS/weight solve and output:4096 cases3.21-3.38s;71080 MPI1 4.432s;
MPI4 sharing one Ada12.333s; perturbed4.604s. This is correctness evidence,
not a distributed speedup benchmark. No long evolution was used for these
tests. Authoritative results: `enriched_v5/{manifest,coupled_analysis}.json`.

## Angular balance and nonstationary periodic exterior

Full physical angular rate from the P1 momentum moments is -.0341763947051569.
Canonical periodic-face angular transport is -.0341763947051567; difference
1.88e-16. This uses both intrinsic and orbital contributions. Wrapped opposite
faces have different physical coordinate representatives, so angular transport
across the periodic cut must be included; it cannot be forced to zero.

The original radial IC is stationary locally in the interior, but not globally
on this periodic square. After subtracting geometric moment transport,
normalized RHS maximum is .4127 in the outer region while r<.047 is below
4.0e-12. Each cell's actual exterior numerical flux is kept. The unchanged
native IC's long-time boundary response must be separated from the final
GIZMO-matched setup/BC test; no equilibrium forcing of the exterior is allowed.

## Thermal tangent gate: current fixed-reference space fails

Independent C++ reference in `enriched_v5_reference_v2` reads the actual native
basis/coefficients and polygons, uses32x16 and64x28 positive volume rules,
and solves the actual Gram matrix. It compares

    native dU/dt coefficients = G^-1 [Mdot - (Gdot)c]

with projection of the independent analytic full Euler strong derivative.
The geometric term is the face integral of w_n U phi, evaluated separately.
This is an instantaneous semidiscrete test, not a trial time step.

For the uniform velocity perturbation b=(.013,-.007) at fixed potential,
the exact initial derivative is -b.grad(T_b U_eq); initial div(v)=0 and P
is constant, hence the physical initial Pdot=0. The common fixed scalar space
contains U_eq, but need not contain its translated tangent. Independent
projection and actual native RHS reproduce the same spurious Pdot. This
separates representation error from the screened solve or missing face terms.

Pressure derivative below is divided by P0, in inverse code-time units:

| Cell | native max abs(Pdot)/P0 | volume L1 abs(Pdot)/P0 | volume mean Pdot/P0 |
|---|---:|---:|---:|
|67943, origin corner|13735.08|703.503|+.0535379|
|68000, origin corner|45783.59|2339.066|-.178460|
|68527, origin corner|45783.59|2339.066|+.178460|
|68584, origin corner|13735.08|703.503|-.0535379|
|67935, old depletion|1.382415|.266094|+3.219e-6|
|68001, old depletion|.880719|.110295|+1.539e-5|
|68526, old depletion|.880728|.110295|-1.539e-5|
|68592, old depletion|1.382405|.266094|-3.219e-6|

The largest point values occur at vertices, but volume quadrature also reaches
45590.84 and volume L1 reaches2339.066. The problem is not confined to a
measure-zero origin. Independent projected volume L1 agrees to about1e-7
absolute or better in these cells; changing quadrature32x16 to64x28 leaves the
large values unchanged. Unprojected exact thermal derivative is <=2.23e-10
over these samples. Native/independent moment RHS discrepancy is <=1.82e-11
of area*component_scale*(speed/covariance_length) in the core.

Without perturbation, selected central native volume L1 pressure derivative
is at most7.45e-7/P0-time and cell means about1.5e-10/P0-time. The severe
perturbed error therefore cannot be dismissed as initial-pressure roundoff.
Initial-state reconstruction and conservation alone are insufficient for the
physical thermal accuracy needed by the active goal. No limiter or extra heat
was applied to conceal this failure, and no evolution from this candidate is
started. Existing code, mass/weight solver and failed records remain available
for the revised operator.

## Next bundled review

Completed review: `ASTRA_NODAL_SBP_PLAN.md` selected positive-diagonal nodal
Euler state and a kinetic/pressure-compatible split volume flux. The actual
implementation now has a separate rejection in `NODAL_SBP_RESULTS.md`:
thermal chain cancellation improved, but cold divergence accuracy and the
advected density-wave local stability gate fail. Preserve both candidates'
evidence; neither is used for evolution.

Historical review request: choose the next implementable thermodynamically compatible
representation/operator. Adding translated tangent functions is a possible
diagnostic, not an automatic solution: derivative enrichment near the literal
origin cusp raises regularity/conditioning questions and cannot replace a
generic Euler/shock consistency argument. A moving representation must include
its true test-function time derivative and cannot mix unrelated bases in SSP.

Kinetic/entropy-compatible split forms are another candidate, but a pressure
or entropy label alone does not establish local linear stability. Relevant
primary literature to assess, without importing their theorems to our meshes:
[Chandrashekar](https://arxiv.org/abs/1209.4994),
[Chan](https://arxiv.org/abs/1708.01243),
[Schnuecke et al.](https://arxiv.org/abs/1812.09093), and
[Ranocha and Gassner](https://arxiv.org/abs/2009.13139).
Their SBP/quadrature/moving-mesh assumptions require an explicit construction
on actual Laguerre polygons. No residual forcing, independent kinetic
reservoir, thermal refill, IC smoothing or invented Re is authorized.
