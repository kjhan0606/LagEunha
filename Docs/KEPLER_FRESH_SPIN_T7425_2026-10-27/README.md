# Affine spin/strain postprocessing at fresh-weighted t≈7.425

CPU job **437596** passed the analytic square, rigid rotation, simple shear
and local Kepler checks, then analyzed the accepted 49,038-cell disk state.
It reconstructs polygon mass moments assuming uniform subcell density and
projects the saved affine velocity gradient into rigid spin plus orthogonal
non-rigid strain. The [result](weighted.json) finds the rigid-spin component
is `11.3293%` of the reconstructed affine kinetic energy. At `t=6.5`, the
same analysis gave `11.3290%`; the change is negligible at the precision of
this reconstruction. The gradient RMS relative to analytic Kepler gradient
is `6.34e-4`.

This result shows that a persistent kinematic rigid-rotation component exists
inside the cells, but the native code has not evolved it as a separate spin
state, torque, or energy reservoir. Subtracting the reconstructed rigid-spin
energy from the current thermal energy while arbitrarily keeping the old
total-energy reference fixed would leave nonpositive thermal energy in about
`15.58%` of disk mass. That is a counterfactual constraint, not a proposed
initialization. A physical spin model still needs its torque and energy-work
law and a consistent initial energy reference.
