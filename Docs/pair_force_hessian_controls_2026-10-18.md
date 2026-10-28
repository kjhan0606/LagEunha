# Force derivative confirms negative thermal curvature; pressure controls

On the same actual t10 geometry and pair 156/157, recomputed Voronoi ridges
and evaluated the original reversible TCF pressure-force formula after
positive/negative center-of-mass-preserving perturbations. Each pressure
model evolves adiabatically with its own fixed cell entropy under displacement.
The restricted stiffness is minus the derivative of the projected pair force.

At epsilon/gap=.0005:

| Baseline pressure model | Minimum symmetric stiffness eigenvalue |
|---|---:|
| Actual evolved pressure | -1.90459517e-4 |
| Uniform pressure | +2.60018988e-5 |
| Initial entropy on evolved volumes | -2.07195837e-6 |
| Actual minus initial-entropy pressure | -1.88387986e-4 |

Actual force derivative differs from the independent energy finite-difference
Hessian by at most 3.49e-10. Halving epsilon reduces antisymmetric residual
from 8.70e-10 to 2.18e-10, consistent with centered finite differencing.
H_actual=H_initial_entropy+H_entropy_excess closes within 2.62e-16. Individual
eigenvalues do not generally add and are not causal percentages.

The uniform-pressure control on the same deformed geometry is positive in
this two-dimensional pair subspace. Small negative curvature is already
present for initial entropy; accumulated entropy-associated pressure changes
strongly increase the observed negative curvature in this state. This is
consistent with pressure/geometry feedback, rather than a force/energy sign
inconsistency in these independent diagnostics.

Limits: one selected pair, restricted displacement subspace, moving state
not an equilibrium, no gravitational Hessian, gyroscopic terms, dissipative
Jacobian or global spectrum. It does not prove an unstable full-system mode,
that entropy originally caused pairing, or that removing heat stabilizes
the actual scheme. No production force or physical trajectory was modified.

Next: establish whether this negative direction persists when gravity and
the instantaneous circular-frame dynamics are included, and whether it is
representative of the problematic pairs. Any cure must retain energy work,
shock heating and the approved weight reduction rather than clip this term.

Evidence: canonical_pair_force_hessian_controls_2026-10-18.json.
Script: audit_pair_force_hessian_controls.py.
Remote results: /gpfs/kjhan/LagForce/Codex_pair_force_hessian_controls_20261018.
