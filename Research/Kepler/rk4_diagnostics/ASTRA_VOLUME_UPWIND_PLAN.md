# Astra next bounded choice: reconstructed Euler HLL volume fluctuation

Completed implementation: `VOLUME_UPWIND_A_RESULTS.md`. The concrete P2 SBP
gradient/midpoint reconstruction fails A (actual angular weak balance) and
also damages the cold pressure tangent. This candidate stops before C/P4 or
runtime integration, as required below. This plan remains as the historical
decision and acceptance contract; it is not a live approval for evolution.

2026-10-21. One bundled evaluation of `NODAL_SBP_RESULTS.md` and recorded
v2/v3/v4 matrices. Both split KEP/PEP and central Euler control are excluded
from evolution. Driver codes/runs. The next candidate is not certified.

## Candidate

For canonical algebraic node pair a,b with fixed skew SBP direction and
measure, use the same ALE metric in all terms:

    f_new = f_KEP(Ua,Ub) + f_HLL(Ua_reconstructed,Ub_reconstructed)
                           - f_KEP(Ua_reconstructed,Ub_reconstructed).

Reconstruct rho,vx,vy,P in the same state-independent polynomial space and
form trace total E from EOS. These are flux traces only; never replace the
stored nodal conserved state, mean E or thermal energy. Constant v/P density
waves and exact affine velocity must retain their trace properties. No
adjustable damping factor. A limiter must not be assumed to repair the smooth
linearized mode when it is inactive. Invalid traces are explicit failures
until a conservative reconstruction limitation is derived.

Driver's concrete first construction: each node uses the same P2 SBP
derivatives of the primitive nodal field. Extrapolate each side linearly
to the common physical midpoint. Exact affine primitive fields then give
matching traces. HLL acts on all four Euler fields, with endpoint and Roe
acoustic speed bounds; neither the bounds nor this high-order formula is
advertised as a universal positivity or entropy theorem.

## Sequential gates: A, then C, then B

A. Existing P2/24 actual polygon patch first: pair antisymmetry and all
Euler conservation, vanishing correction for affine matching traces,
boost/ALE consistency, and actual angular weak balance. Use

    Lh=sum_a M_a x_a cross (rho v)_a

including moving-node position derivatives; compare physical boundary angular
flux and gravity torque. Report volume correction angular contribution
separately. Pair momentum conservation alone is insufficient. Do not create
spin to absorb an alias residual or erase it with a torque correction.
If A fails, STOP this candidate before C, P4 or runtime integration.

C. If A passes, same advected w=0 density wave as v3, actual analytic JVP
including reconstruction and full numerical finite differences. A merely
smaller positive eigenvalue is not a pass; retain limiter branch meanings.
The favorable w=v subcase is insufficient.

B. Only if A/C pass, one P4 positive nodal SBP construction, not a degree
sweep. Same24 cells and b=(.013,-.007). Separately report KEP Pdot, upwind
correction Pdot, -gamma P div_h(v), total residual and original nonpolynomial
IC mass-quadrature discrepancy. P4 is aimed at divergence accuracy, not an
assumption that Riemann dissipation corrects it. If heating spoils cold
accuracy, stop without tuning the correction's amplitude.

After all three only: nonaffine Laguerre GCL, screened weight solve, native
CUDA/MPI runtime/restart, conservative positivity, Sod/Sedov, and refined
Kepler .2 then t6-10 then longest GIZMO duration. No early long evolution.

Euler wave-speed dissipation is numerical, not imposed finite physical Re,
but its actual cold heating must be measured. Do not transfer graph-viscosity
or convex-limiting theorems to this corrected high-order operator without
their assumptions. [Guermond et al.](https://arxiv.org/abs/1710.00417).
