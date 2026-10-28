# Restricted rotating-pair tangent and gravity controls

For t10 pair 156/157 use reduced mass mu=4.8806903843e-6 and a frozen frame
rotating at the circular frequency of its mass center, Omega=6.9026180574.
Project the central-potential Hessians onto the pair relative coordinates,
add thermal stiffness, subtract mu Omega^2 I, and include Coriolis 2 Omega J:

    qddot + 2 Omega J qdot + (Hthermal+Hgravity-mu Omega^2 I)q/mu=0.

This constrains pair mass center and all surrounding generators. It freezes
the evolved geometry and pressure, omits dissipation and entropy perturbations,
and is not the full moving-trajectory variational equation.

Maximum eigenvalue real parts in this restricted model:
actual pressure 5.7608, uniform pressure 2.1395, initial-entropy pressure 2.6128.
These are not measured simulation growth rates or evidence of a global mode.
In particular the positive uniform-pressure control prevents attributing this
model's growth solely to the actual state's negative thermal-energy curvature.
Rotation/tides and constrained pressure response must be treated together.

Pressure-free finite-pair control has purely imaginary eigenvalues (to about
1e-15 real parts), with frequencies 11.7066 and .02890. The coincident local
circular control gives two neutral eigenvalues and frequency 11.7065, matching
kappa^2=4 Omega^2+r d(Omega^2)/dr=137.04217>0. Normalized eigenvector residuals
are below 3.1e-16. Thus no unconditional gravity-only exponential growth is
introduced by these basic controls; this still does not validate the constrained
pair model as a complete instability explanation.

Next required gate: release surrounding-particle and center-of-mass degrees
of freedom, or directly measure paired perturbed trajectories using the actual
solver. Include damping and evolving geometry before proposing a correction
based on these eigenvalues. Do not equate positive stiffness alone with
stability in the shearing rotating problem.

Evidence: canonical_pair_rotating_tangent_2026-10-18.json and
canonical_pair_gravity_controls_2026-10-18.json.
Script: audit_pair_rotating_tangent.py. No physical trajectory changed.
