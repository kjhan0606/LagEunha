# QBP signal decomposition on actual evolved t10 geometry

The prior QR operator export was from a circular-velocity diagnostic at t2
(maximum departure from local circular speed 4.22e-15). It cannot represent
evolved t10 perturbations. A new terminal geometry/force export therefore
read the preserved unmodified t10 checkpoint 2305. No evolution or steering
was applied. Original long trajectories remain untouched.

The actual C QR gradient map reproduces stored gradients within 1.55e-13.
All fitted statuses are valid and arc flags zero. Reconstructed QBP delta and
qin match printed source outputs at relative L2 1.57e-7 and 1.35e-7.
At fixed geometry, v=v_circular+v_perturbation gives additive delta within
8.89e-15; compression switch chi is nonlinear and is recomputed separately.

| Current-radius region | Actual qin RMS | Circular qin RMS | Actual mean chi | Circular mean chi |
|---|---:|---:|---:|---:|
| inner r<.5 | .025346 | .000332865 | .13617 | .0006687 |
| edge .5<=r<.6 | .00934183 | .000092636 | .04032 | .0002227 |
| core .6<=r<1.9 | .0000565934 | .00000652011 | .00009253 | .000008136 |

Statistics are face-area weighted, directed records, current-radius bins.
Inner delta RMS is .044580 actual, .0003032 circular, .044623 perturbation.
The actual t10 compression signal is thus dominated by evolved velocity
perturbations, not the residual response to circular velocity on that same
geometry. This does not identify the initial perturbation source, prove that
TCF is necessary/sufficient, or define an additive heat fraction. No minmod,
FACE_Z heating or evolved counterfactual was used for this conclusion.

This strengthens the reason to address mode generation upstream rather than
simply reduce late-time AV. A circular-background-subtracted viscosity is not
adopted: it would need independent shock/perturbation consistency checks and
must not steer velocities. The geometric versus material volume mismatch
remains a candidate source requiring causal early-time testing.

Evidence: canonical_qbp_t10_decomposition_2026-10-18.json.
Script: analyze_qbp_t10_decomposition.py.
Remote raw operator, cache, face/RHS files and exit code:
/gpfs/kjhan/LagForce/Codex_qbp_t10_decomposition_20261018.
