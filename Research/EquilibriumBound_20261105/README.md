# Minimum uniform-pressure equilibrium stiffness control

See ../../Docs/UNIFIED_HYDRO_EQUILIBRIUM_RETIREMENT_2026-11-05.md.
The same native cells and fixed masses define a rest state with uniform
pressure. Constant pressure annihilates the geometric force derivative, so
the positive acoustic block is the complete inviscid pressure linearization
of the proposed semidiscrete model at that equilibrium.

B uses analytical rectangular-matrix norm inequalities rather than treating
a Ritz value as a maximum-eigenvalue upper certificate. C reuses the previous
independent Rayleigh value and pressure-range comparison. With explicit
floating padding the frequency-ratio lower estimate remains 36.02. This is
not an interval-arithmetic certificate or a whole-solver/GIZMO runtime claim.
No new eigensolve, coefficient solve, or PDE step was performed.

The existing local psi_i times fitted-strain construction is retired from
production consideration under the bundled Astra stop rule. Source, proof,
failed attempts, and the original native trajectories remain preserved.
The result does not reject all C=B+LR models and does not complete the full
physical-accuracy, dimensional-suite, or matched-speed objective.
