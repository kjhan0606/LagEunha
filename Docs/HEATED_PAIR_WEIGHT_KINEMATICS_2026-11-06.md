# Heated material pair: pressure-driven boundary, geometric weight rate

User example: initially equal material cells, stationary generators and
midpoint face; inject heat into the left cell. Reused the existing native
header Exam/gfs_exact_riemann_ale.h, not a new Riemann implementation.
Artifact: Research/MaterialVolume_20261105/heated_pair_material_20261106_v1.
audit.c compiles as C11 and exits0 for seven heat fractions in1D and the
same normal problem embedded at n=(.6,.8) in2D. No3D test or PDE run.

Normalize rho=1, initial P=1, gamma=5/3 and generator separation d=1.
At fixed initial volume the internal-energy increase fraction equals the
pressure increase fraction. The existing exact Riemann solution supplies
contact pressure p* and speed u*. Evaluate its ALE flux at w=u*:

    mass flux = 0,
    momentum flux = p* n,
    total energy flux = p* u*.

All checked to floating precision (largest absolute energy-work discrepancy
5.69e-14 for the strongest heating). Equal states give exactly zero motion.
For10% heating p*=1.04976187664, u*=.0378003207256.
At W=0 and stationary generators, Delta Wdot=2d u*=.0756006414512.
The2D rotated normal problem gives the same scalar p*,u*.

## Geometry, including moving generators

In1D, face b=(x_i+x_j)/2+(W_i-W_j)/(2d), d=x_j-x_i.
To represent a prescribed material contact speed u* exactly for this pair:

    Delta Wdot = 2d [u*-(v_i+v_j)/2] + (Delta W/d)(v_j-v_i).

Only weight differences matter. In multiple dimensions at face centroid,
let n=(x_j-x_i)/d and c=x_face-(x_i+x_j)/2. The corresponding pair requirement:

    Delta Wdot = 2d [u*_n-(v_i+v_j) dot n/2] + 2(v_j-v_i) dot c.

This is a kinematic conversion, not an independent weight equation of motion
or a need for auxiliary kinetic energy. Geometry and physical Riemann work
must use the same actual boundary velocity.

## Limits that must not be hidden

The earlier two-cell adiabatic pressure-balance volume ratio
[(U0+Delta U)/U0]^(1/gamma) is conditional on both cells preserving entropy.
It is NOT a general prediction of final static equilibrium after heating a
finite domain. The actual discontinuous heating problem creates a rarefaction
on the hot side and a shock on the cold side. For10x pressure, right-star
P/rho^gamma rises by22.74%; the per-cell isentropic final-volume formula does
not apply to that shocked material. Pressure equality at a moving contact is
also not global static equilibrium. Finite-domain reflections/energy
partition require time evolution and boundary conditions.

For many cells, independently computed edge Delta Wdot need not sum to zero
around graph cycles. A scalar Wdot per cell therefore cannot generally realize
all independent Riemann contact speeds while prescribed generators move.
Projecting only net volume does not erase this mismatch. A multidimensional
compatible boundary/pressure solve must explicitly handle it, conserve work,
respect the geometric admissibility bound, and meet runtime constraints.
Do not simply impose pair rates, clip weights, or discard the work defect.

The existing1D solver remains the physical reference. This diagnostic supplies
a concrete pressure->material speed->weight-rate link without extra inertia;
it does not complete its2D implementation or qualify shock behavior of the
smooth native variational candidate. All prior trajectories remain preserved.
