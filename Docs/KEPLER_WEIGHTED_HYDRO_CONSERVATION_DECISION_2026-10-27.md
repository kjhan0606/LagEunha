# Weighted Kepler hydro: conservation decision after actual t=6.5 probes

## Measured constraint at the live solver state

The accepted bounded q2 aligned and 5%-jitter native CUDA+MPI+OpenMP
checkpoints at exact `t=6.5` both pass 17 physical gates. Static Laguerre
weights recover each of the 49,038 disk cells' initial areas to maximum
relative error below `3.3e-11`. A 16-iteration weight-rate solve is close to
the converged volume-rate solve for one current `dt≈0.00194` geometry step.
These are genuine evolved native generator positions, with no site steering.

The [material-rate projection](KEPLER_Q2_BOUNDED_MATERIAL_RATE_T65_2026-10-27/README.md)
then tests whether weight rates can also make every disk face move with the
analytic Kepler fluid. They cannot on these two snapshots: the least-squares
optimum still retains 86.6% aligned and 88.4% jittered of the fixed-weight
RMS face-slip magnitude. Compared with the almost volume-preserving solution,
it lowers slip RMS by only 1.49% and 1.14% while restoring 42.9% and 29.6%
of the original cell-volume-rate L2. One scalar rate per cell is therefore
not a useful standalone route to material faces at these states. The
projection uses analytic face velocity, not native Riemann reconstruction.

The later [t=10 native comparison](KEPLER_Q2_BOUNDED_T10_DIAG_2026-10-27/README.md)
adds a sharper problem target: aligned rings pass all 17 checks, whereas
5%-jitter rings fail the `disk_ballistic_residual < 1e-3` gate at
`1.19668e-3`. Jittered disk pressure-acceleration RMS rises monotonically
from `2.51e-7` at t=6.5 to `3.88e-7` at t=10, predominantly tangential;
aligned force stays near `1.7e-8`. The reconstructed force agrees with the
existing GPU 64-cell dump to at most `1.01e-11` in acceleration vector.
A deliberately unphysical one-step pressure-force-off CUDA control from
the accepted jitter t=6.5 state changes velocity by `2.50690e-7 dt`,
matching the independently reconstructed `2.50685e-7 dt` pressure force.
Its [t=10 continuation](KEPLER_Q2_PRESSURE_FORCE_CAUSAL_2026-10-27/README.md)
tests whether this force explains the accumulated residual; its result
must not be called a physical hydro solution.

Static target-volume weights still exist at the **failed** jitter t=10
positions: maximum relative area error `1.78e-10`, maximum
`|delta W|/d²=0.0171`, no faces beyond the current `0.8` guard. A frozen
counterfactual that changes only weighted face geometry leaves pressure
force nearly unchanged; also recomputing pressure adiabatically from the
weighted cell areas lowers jittered disk pressure acceleration RMS
`3.88e-7→1.05e-8`. This identifies the area/EOS pressure variation as the
main *instantaneous* force source at t=10 and makes stage-coupled weighted
geometry a concrete candidate. It does not solve the dynamic energy and
mass-transport questions below.

## Exact pressure-work identity

Write a cell's power area as `V_i(x,W)`. For a finite face `i|j`, let
`g_ij=L_ij/(2d_ij)`. Its weight-driven area rate is

`Vdot_i|W = Σ_j g_ij (Wdot_i-Wdot_j)`.

With fixed mass and adiabatic internal energy, `Udot_i=-P_i Vdot_i`.
The current variational site-pressure force is the derivative with respect
to **site positions at fixed weights**. Its mechanical power is
`Σ_i P_i Vdot_i|x`. Therefore, if that force is retained while weights move,

`d(K+U)/dt = -Σ_i P_i Vdot_i|W
             = -Σ_(i|j) g_ij (P_i-P_j)(Wdot_i-Wdot_j)`

apart from gravity and boundary terms treated consistently. This term
vanishes for frozen weights or exactly uniform pressure, but it is not an
identity for a general weighted compressible flow. In the actual cold
`t=6.5` states it is small: the native-pressure signed right-hand face sum
is `-1.240e-13` aligned and `+6.310e-12` jittered per unit time; using an
adiabatically reweighted active-disk pressure proxy gives approximately
`-6.37e-13` and `-6.22e-13`. Face and cell sums agree to roundoff in the
[recorded audit](KEPLER_Q2_BOUNDED_MATERIAL_RATE_T65_2026-10-27/README.md).
Smallness in this near-uniform-pressure test is not a general energy law.

## Why the approved W=0 reduction matters

The approved reduction is: zero weights give ordinary Voronoi geometry and
fixed cell masses. A full arbitrary-moving-face ALE flux would instead have
face mass rate proportional to `rho_face L (v_fluid,n-u_face,n)`. Actual W=0
Voronoi faces in Kepler shear have nonzero normal slip, so a full ALE flux
would change masses even at W=0. Thus a single unmodified ALE flux cannot
both be exact on these faces and reproduce the approved W=0 fixed-mass
method. This is a discrete-model choice, not a GPU implementation detail.

At the weighted `t=6.5` geometry, a frozen-state flux proxy over disk-disk
faces gives gross absolute crossing of `0.010997` aligned and `0.015203`
jittered disk masses per physical time unit. Pairwise mass rates cancel
globally; the p99 **net** cell change over one current step is only
`6.92e-7` and `6.13e-6` of cell mass. Gross crossing can still exchange
tangential momentum and energy across the shear. The proxy omits disk
boundaries and native reconstructed states, and it cannot be extrapolated
to `t=600`.

For the same proxy, transporting `mdot v_K(c)` at the face centroid gives
an absolute sum of off-center orbital face torques `6.032e-4` aligned and
`1.642e-3` jittered per unit time. The signed sums nearly cancel globally,
but a cell-spin ledger is required for **pairwise** angular-momentum
accounting if such fluxes are used. The pair identity is
`(x_j-x_i)×F_p -(c-x_i)×F_p +(c-x_j)×F_p = 0`; its measured residual is
below `5.3e-23` in the [proxy](KEPLER_Q2_BOUNDED_MATERIAL_RATE_T65_2026-10-27/README.md).
This is an accounting test, not evidence that physical cell rotation
removes shear heating.

## Native implementation gate

1. Preserve the current `W=0`, fixed-mass, variational RK4 solver byte for
   byte as the reference branch. Measure both reference and candidate from
   the **same accepted checkpoint**; never restart from the original IC to
   compare late times.
2. Attach a device-side weight-rate solve to the already resident FaceCSR;
   measure 8, 16 and 32 iterations and the true whole-RK4-step cost. The
   target volume rate must come from a physical divergence prescription,
   with zero as the Kepler equilibrium special case. Enforce positive power
   cells and the native candidate/containment guards at every RK stage.
3. Before activating weighted hydro, select a face update with an explicit
   ledger for mass, orbital momentum, internal energy, gravitational work,
   and any energy stored in the weight degree of freedom. A standard ALE
   update is a conservative comparison branch, but its W=0 mass behavior
   differs from the approved reduction. A fixed-mass weighted update needs
   an explicit correction for the pressure-weight work above; merely
   inserting `Wdot` into `PdV` breaks the existing site-force energy
   identity. Any weight-linked mass exchange must use antisymmetric face
   fluxes and demonstrate why it improves the actual density/entropy error.
4. Start with a one-step native `t=6.5` restart, same sites and pressure
   field. Compare W=0 reduction, volume/face-speed/work identities,
   pairwise conservation, entropy sign, density noise, angular momentum,
   and GPU wall time. Then cross `t=10` and extend toward the GIZMO long
   endpoint only if these gates pass. Spin is a separate physical state:
   its off-center torque and rotational energy need a pairwise ledger;
   physical Kepler viscosity remains `Re=∞`.

The present evidence supports **geometry feasibility** and rejects
weight-rate-only exact face materiality at `t=6.5`. It does not yet select a
conservative weighted hydro operator or prove long-time Kepler stability.
