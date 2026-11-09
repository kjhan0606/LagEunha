# Acoustic residuals with matched adjoints: algebraic screen, not adoption

## Why this screen

The old coupled-weight acoustic screen converted geometric face/fluid-speed
mismatch into heat, including in resolved shear. Simply moving to a symmetric
reconstructed central flux retained a checkerboard null. This controlled
comparison separates numerical geometric realization from acoustic residuals,
then pairs each reconstructed residual with its transpose in force/volume work.
It adds no physical state or tunable sensor. It does NOT yet supply a qualified
nonlinear shock solver or native time integrator.

Saved32x32 Gresho geometries (initial, deformed frozen/dynamic) and the SAME
inverse-distance-squared reconstructions as the prior coupled acoustic audit
are reused for a controlled comparison. They are not the face-angular fit of
the more recent KH reconstruction screen; do not mix these results.

## Equations tested

Ti,Tj map cell velocities to left/right normal face traces. Si,Sj are the
corresponding scalar pressure trace maps. Define

    Jv=Tj-Ti, Jp=Si-Sj,
    Zi=sqrt(gamma rho_i P_i), Zj=sqrt(gamma rho_j P_j),
    Z=Zi+Zj, Zeff=Zi Zj/Z,
    Tbar=(Zi Ti+Zj Tj)/Z, C=D Tbar,
    wP=Jp P/Z.

D is signed face-area incidence; Af below is face area. Target volume rate:

    Vtarget = C v + Jp^T Af wP.

Realize it geometrically using ONE graph solve, L Wdot=Vtarget-Bv. Use

    F = C^T P - Jv^T Af Zeff Jv v,
    Udot_i = -P_i Vdot_i + Hi.

For each face assign nonnegative heat to its two cells:

    Hi += Af Zi [wP + (Zj/Z) Jv v]^2,
    Hj += Af Zj [wP - (Zi/Z) Jv v]^2.

Cross terms cancel in the sum. Total heat equals
sum Af [Z wP^2 + Zeff (Jv v)^2], and v^T F+sum Udot=0 because each residual
and adjoint use the same stencil. For fixed mass ideal-gas thermodynamics,
Udot_i+P_i Vdot_i=Hi>=0 supplies nonnegative cell entropy production.

This is an instantaneous semidiscrete identity. Neither positivity during a
finite step nor full nonlinear accuracy follows from it. The high-order
pressure adjoint distributes volume changes through a larger stencil;
Vtarget is generally NOT D(Tbar v+wP). Thus the target is a cell-volume law,
not proof that each geometric face realizes a local material Riemann speed.
That distinction is recorded numerically and cannot be hidden by claiming
zero face mass exchange. No mass flux is evolved in this diagnostic.

## Controlled comparisons

An intermediate split construction used unreconstructed pressure differences
and a separate impedance-weighted global pressure projection. It reduced
geometry-mismatch heating but required two graph solves. Its results are
preserved in split_material_acoustic_20261106_v1. The final residual-adjoint
screen reconstructs pressure too and needs one geometric solve; it does not
claim a native timing result for that solve on these RHS.

The final screen has energy-rate residual below7.22e-16, target-volume residual
below2.23e-16, and nonnegative cell heat in every tested case. Uniform rest is
preserved to roundoff. The regular-grid alternating pressure amplitude has
instantaneous logarithmic decay rate -165.246, so it is no longer invisible.

For cold smooth shear P=1e-6, vx=sin(2pi y), the original geometric-mismatch
model heated at approximately0.468 and0.516 times initial U per time on the
two deformed geometries. The residual-adjoint result is0.035725 and0.028938.
This is an improvement, but still2.9--3.6% of initial U per time; do NOT call
the shear gate passed or a long-time simulation solved. The remaining heat
comes from reconstructed velocity jumps, not pressure residuals. On the
regular grid the same shear heat is roundoff-level.

For initial Gresho, pressure-channel heat drops to3.23e-5 (velocity-channel
heat0.003624); total heat/initialU per time is4.28e-4. The distorted saved
states still have larger residuals. Those states contain preexisting errors,
so this is an operator response comparison, not a trajectory improvement.

## Relation to the verified1D method

Without reconstruction and geometric constraints, the normal-face relations
reduce to the acoustic impedance Riemann formulas. A native C unit comparison
calls the existing exact Riemann solver on a2D face with normal(.6,.8).
Halving small state perturbations quarters the star-speed/pressure errors:
at amplitude.01 errors are7.00e-6 and9.08e-6; at.00125 they are1.09e-7 and
1.42e-7. The two-sided acoustic heat identity also passes.

This verifies the weak-wave limit ONLY. At amplitude1 the differences are
0.0677 in speed and0.1085 in pressure. Therefore this is NOT the exact
nonlinear Riemann method used in the successful1D material tests. A common
strong-shock extension remains required; no1D result was replaced.

## Status

Useful progress is a matched force/volume/heat formulation with an acoustic
checkerboard response and a single geometric solve, plus a quantified remaining
shear error. It is not adopted. Nonlinear shocks, contact/moving-contact time
evolution, weight admissibility, full2D accuracy, MPI/CUDA implementation and
end-to-end cost remain unverified. No3D run started or GIZMO speed pass claimed.

Artifacts: Research/MaterialVolume_20261105/residual_adjoint_acoustic_20261106_v1/:
audit.py/results.json, face_limit.c, executable and face_limit.csv. Python
performs sparse algebra on saved geometry; the face-limit comparison uses
the existing C exact solver. Native production physics was not modified.
