# Actual full-solver response to a pair velocity seed

Four isolated checkpoint copies at t10.0020667693911 received signed relative
velocity perturbations +/-1e-6 and +/-5e-7 on pair 156/157, oriented along the
negative thermal-Hessian eigenvector. Mass-weighted coefficients preserve
pair total momentum (measured residual below 9e-22). Positions, mass and
internal energy were not changed. This is an initial diagnostic perturbation,
not steering during evolution. Every other particle evolves normally.

All four runs completed eight actual HG/RK4/HG steps, dt=.0005, to
t10.006066769391104 with finite accepted states and zero retries. A previous
same-binary unperturbed branch supplies the central reference. The central
velocity derivative's mass-weighted norm gain is .90406227 for both amplitudes:
about 9.6% reduction over this interval. Half/full derivative difference is
7.16e-8 relative; even departure from reference halves from 9.63e-6 to 4.82e-6
relative to odd response. About .196% of squared mass-weighted velocity
response lies outside the selected pair at this endpoint.

This rules out immediate velocity amplification for this particular seed
over this interval, despite negative thermal curvature in the restricted
model. Full forces, dissipation, neighbor motion and evolving geometry matter.
It does not rule out another growing seed or later growth: this was a velocity
seed, not a complete unstable eigenvector in position/velocity/entropy space.
Velocity norm is not total perturbation energy, and .004 is not a long-time
physical validation. The restricted eigenvalue 5.76 must not be presented as
the measured growth rate of the simulation.

Next: continue these same diagnostic branches from their final checkpoints
to a longer common interval if needed, rather than recreating them, and
include pressure/position response in the norm before inferring stability.
Do not alter the original long trajectories or remove dissipation based on
the pressure-only Hessian.

Evidence: canonical_pair_tangent_branches_2026-10-18.json.
Scripts: run_pair_tangent_branches.py, analyze_pair_tangent_branches.py.
Remote data: /gpfs/kjhan/LagForce/Codex_pair_tangent_branches_20261018.
