# Positive-diagonal nodal SBP: identities improve, two physics gates fail

Subsequent evidence and rejected next candidate:
`VOLUME_UPWIND_A_RESULTS.md`. That report also measures the original split
flux's angular weak residual; it is not certified by the SBP identities here.

2026-10-21. Driver implements/tests; one completed spatial bundle for Astra.
Authoritative source: `/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`,
origin `git@github.com:kjhan0606/LagEunha.git`, HEAD `8e25dfe0` plus retained
working changes. New native helpers: `Exam/gfs_sbp_flux.h`,
`Exam/gfs_sbp_polygon.h`. No legacy runtime flag or timestep changed.

## Construction and scope

Same actual original failure polygons and canonical-face direct neighbors:
24 cells, 36 nodes each. Positive diagonal M, boundary nodes selected directly
from the volume set. Centroid-fan cubature uses Gauss3 along the face and
3-point Radau for weight t dt (including t=1). Its volume polynomial degree
and face degree satisfy the P2 SBP compatibility contract. No vertices sample
the narrow original origin density cusp with a finite lumped endpoint weight.
No change to the continuous rho,v,P, gravity, epsilon or original outputs.

Actual state is 4 conserved Euler values per node, with no enriched-state
projection or dense state compression. For polynomial Vandermonde V and
P=(V^T V)^-1 V^T, A=M V_k-E_k V/2, construct

    S=A P-(A P)^T+P^T (A^T V) P, Q=S+E/2.

Only the measured floating-point skew part of A^T V is used; unmodified
cubature compatibility is exported. This constructs an operator independent
of the fluid state, not an equilibrium subtraction. Pairwise skew volume
fluxes are opposite and the boundary uses the full existing Euler HLLC/HLL
flux. Boundary and volume use the same nodes. Same-node gravity has zero
thermal source. The current matrix supports an affine mesh-motion GCL;
actual nonaffine Laguerre vertex kinematics and runtime MPI/weights remain
unconnected. This is explicitly not a complete moving-cell gate A pass.

Frozen evidence: `nodal_sbp_v1/` (finite difference only), `nodal_sbp_v2/`
(stationary translating wave with independent analytic derivative), and
`nodal_sbp_v3/` (advected wave, generic native analytic volume JVP plus exact
matching-trace HLLC variation). C++ performs all flux/RHS/Jacobian work on
user-authorized lageunha; Python only diagonalizes recorded matrices and
summarizes outputs. Runs are terminal and no evolved state exists.
`run_nodal_sbp_v1.sh VERSION [advected]` refuses an existing output directory.
Sources/binaries and SHA checks are frozen separately for every version.

## A: algebra and actual geometry

12000 two-point samples: symmetry 0, Euler ALE consistency 4.57e-15,
kinetic identity 1.87e-14, boost 6.61e-15 (normalized absolute metrics).
Positive radial rule moments through degree4: error 5.55e-17.

Across the24 actual polygons:

| Quantity | Maximum error |
|---|---:|
| Q+Q^T-E, normalized | 0 |
| Q*1, normalized | 8.27e-17 |
| P2 derivative exactness, normalized | 1.26e-16 |
| unmodified cubature compatibility, normalized | 5.93e-16 |
| affine metric div_h(w)-trace(B) | 2.02e-14 |
| same-node gravity thermal contribution / P0 | 0 |

Minimum positive node mass / cell area .0115741. Initial pressure relative
error <=5.65e-13. Cubature is not exact for the nonpolynomial original IC:
max initial cell-mass discrepancy against prior integrated native IC is
.00121123 relative (origin-adjacent cells); max total-energy discrepancy is
.000362777 of that cell's thermal energy. These are disclosed initialization
errors, not quantities repaired by changing the heat or density.

## B: thermal chain identity passes; Kepler divergence accuracy does not

Original equilibrium, perturbation b=(.013,-.007), and the same perturbation
with a prescribed affine mesh velocity are evaluated separately. For exact
constant initial P and div(v)=0, physical initial Pdot=0. Numerically the new
operator gives Pdot=-gamma P div_h(v), with max residual/P0=1.51e-8/time.
This eliminates the earlier dense-projection kinetic/thermal chain defect
in this spatial check, without a heat reservoir or limiter.

However max abs(div_h(v))=2.23874649/time and max abs(Pdot)/P0=3.73124487/time.
Four origin-adjacent cells have volume L1 abs(Pdot)/P0=.339445/time despite
near-zero signed mean. Neighbor cells reach L1=.383060/time. Cancellation
of a cell mean is insufficient. This remains a failed cold Kepler accuracy
gate, not a success inferred from the improved identity. Do not run it long
or force div(v)=0 by changing the physical state.

## C: advected density wave rejects this volume flux/operator candidate

Periodic 2x2 square polygons, same P2 construction, 144 nodes/576 conserved
variables. rho=1+.6 sin(2*pi*x), v=(.8,-.3), P=1, gamma=1.4. Each of24 face
nodes is evaluated once canonically with opposite signs. Max summed
conservative rates <=4.10e-16. No time evolution or residual subtraction.

1. With prescribed uniform translating metrics w=v, the density wave is a
   stationary Euler ALE state. Its exact primitive acoustic Jacobian has
   max real eigenvalue 1.21e-14 and max energy-weighted symmetric eigenvalue
   5.20e-14. Entropy-variable columns are identically zero. Full-state finite
   differences show spurious near-zero/Jordan eigenvalues shrinking with the
   differentiation interval; these are not called physical instability.
2. This favorable subcase is insufficient. With w=0, the wave advects through
   the same nodal operator. The exact conservative Jacobian has **max real
   eigenvalue 1.81039321666/time**, min real -81.2889530, max abs imaginary
   51.3332588. Native finite differences confirm it:

| relative perturbation | max real eigenvalue | max matrix error vs exact |
|---|---:|---:|
| 1e-4 | 1.81039284578 | .00301939 |
| 1e-5 | 1.81039318197 | .000301937 |
| 1e-6 | 1.81039321353 | .0000301986 |

The analytic JVP differentiates the actual shared volume flux. At identical
base traces the HLLC derivative uses the exact linear acoustic/contact fan,
with the actual mesh velocity fixed. This matches the full numerical flux
finite differences. The advected base has a nonzero mass RHS; no claim is
made that its complete nonlinear long-time solution is determined by the
frozen Jacobian alone. Nevertheless this is exactly the local stability
screen requested by Astra and the known issue in the primary literature.
It is not cured by an HLLC boundary or by KEP/PEP labels. The stationary
translating wave is not substituted for the actual advected check.

[Ranocha & Gassner](https://arxiv.org/pdf/2009.13139), Sections3.3-3.4,
also show that the Shima arithmetic-density KEP/PEP flux can admit locally
growing perturbations despite an unperturbed density wave running a long
time; surface HLL dissipation alone did not eliminate their example's
growth. Those paper meshes/results are supporting context, not a theorem
for these particular polygon operators. Our rejection uses our matrices.

## Next decision / retained final target

### Diagnostic control isolates the smooth volume-flux issue

`nodal_sbp_v4/` keeps identical polygons, M/Q/S/R, boundary nodes and the full
HLLC boundary flux. Only its smooth volume flux is replaced by the arithmetic
average of complete physical Euler ALE fluxes. This is a diagnostic control,
not a selected candidate. Its exact volume derivative and matching-trace
HLLC derivative are again checked against the full numerical flux.

The same advected density-wave Jacobian now has max real eigenvalue
1.23346e-14; native finite differences at1e-4,1e-5,1e-6 give max real values
3.47e-11,4.81e-10,4.63e-9 respectively. This identifies the positive1.81039
mode as dependent on the chosen split volume flux, not an unavoidable defect
of these polygon SBP matrices. It is consistent with the arithmetic full
Euler control reported in the primary paper.

The control fails the other required property badly: on the same24 original
Kepler cells max thermal-chain residual/P0 is20977.53/time, and max absolute
Pdot/P0 is20973.83/time. Stable density-wave spectra alone cannot select a
cold-flow operator. Neither control nor the split-flux candidate advances
to long evolution. All v1-v4 records and native terminal status are retained.

This candidate is excluded from nodal runtime integration and shock/long
evolution. A bundled Astra review must select the next implementable generic
Euler operator, addressing both local linear stability and cold divergence
accuracy while preserving total energy, momentum, actual angular transport,
same-mesh GCL and shared-face conservation. No arbitrary damping, pressure
reset, exact-residual subtraction, independent kinetic reservoir or generator
steering is introduced. Retain all reusable helpers and failed evidence.
The actual t6-10 interval and longest GIZMO comparison remain the goal;
this spatial bundle proves no long-time Kepler stability.
