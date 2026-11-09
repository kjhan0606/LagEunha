# A physical kinetic metric also requires matching compression energy

Algebraic1D uniform-rest reduction, not a new1D benchmark trajectory and not
a replacement for the validated material Riemann/PLM code. This tests the
proposed polynomial kinetic representation before any native port. Units
h=rho=c=1. For a Fourier mode theta and generalized velocity z=(v,Wdot),
face velocities imply

    D z = wR-wL,       D=(i sin(theta), 1-cos(theta)),
    S z = wR+wL-2v,    S=(cos(theta)-1, -i sin(theta)).

The exact quadratic-profile kinetic integral previously derived gives

    M = diag(1,0) + D^H D/12 + S^H S/20.

The original cell-volume-only compressive stiffness is A0=D^H D. Generalized
eigenvalues of (A0,M) retain a zero branch but the nonzero branch has
omega/theta -> sqrt(6), approximately2.44949, not1. This is independently
verified by the rank-one expression D M^-1 D^H. Thus changing only inertia
would corrupt long-wave sound speed; the energy integral alone is insufficient.
It is a failure of this mixed closure, not a proof that physical internal
kinetic energy is wrong.

## Matched compression representation

The same quadratic displacement profile has strain a+2b*x, with a=D z and
b=3 S z in these units. Its actual linearized compression-energy integral is

    .5 integral (partial_x displacement)^2 dx
       = .5 [|D z|^2 + 3 |S z|^2].

Therefore the matched stiffness is A1=D^H D+3 S^H S. This expression uses
displacement amplitudes rather than velocity amplitudes in the potential;
the same coordinate maps D,S apply. It follows from continuum linear acoustic
energy and is not an empirically tuned restoring term. The compression
integral was checked in1000 random real profiles, relative error <9.05e-16.

The low branch of (A1,M) has sound-speed ratio1.000000000001 at1024-cell
wavelength and1.00000006448 at64-cell wavelength. The missing pressure
variation matters: the original volume energy retains only average strain,
while the physical velocity metric permits nonuniform strain within a cell.
Both fields must use compatible representations.

## Cost and scope

The matched representation introduces a high internal acoustic branch with
omega*h/c -> sqrt(60), about7.746. This is a real extra mode of the chosen
polynomial model, not a measured PDE runtime. It would impose a tighter
explicit harmonic stability limit; do not adopt this as a simple inertia
replacement or claim the <=2x GIZMO goal is met. No independent state was
added to production, no3D calculation was run, and no costly integration was
put into a timestep.

The useful design requirement is now specific: couple velocity and density/
compression representation consistently, and assess whether the fast internal
response can be represented economically without arbitrary damping or a
wrong long-wave limit. Existing affine moments can supply parts of the
kinetic calculation, but do not by themselves supply this pressure closure.
The successful W0 Kepler and1D reference codes remain preserved.

Artifact: Research/MaterialVolume_20261105/weight_metric_acoustic_20261106_v1.
audit.py/results.json retain the failed inertia-only screen;
consistent_pressure.py/consistent_pressure_results.json retain the matched
linear representation and independent integral verification.
