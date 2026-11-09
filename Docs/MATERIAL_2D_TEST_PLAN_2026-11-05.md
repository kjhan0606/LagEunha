# Two-dimensional material-cell validation stage

User approved autonomous2D validation until further instruction, prioritizing transfer of useful1D results over perfecting residual1D peaks. Kepler orbital dynamics excluded from this bundle. RT remains included and uses external gravity. Preserve original production OpenMP/MPI/CUDA, initial conditions, outputs and restart history. No /scratch.

Current state: fixed-mass shared-interface1D PLM/SSPRK2 passed all nine original endpoints, with residual local errors. Existing2D affine rotation/strain implementation is retained. New2D evolution under the same material constraints has NOT been run or qualified. Old2D results must not be presented as new-candidate passes.

## Sequence

1. Uniform flow, stationary/moving density contacts, planar and oblique shocks, smooth acoustics. Verify dimensional reduction and mesh-direction dependence. Planar strip extrusion is only a reduction check, not a genuine2D validation.
2. Cold smooth shear and Gresho vortex. Measure pressure/thermal growth, velocity/density error, angular momentum and resolution dependence. No central Kepler gravity.
3. KH and RT using retained original definitions and physical endpoints. RT includes consistent external gravitational work and initial hydrostatic balance. Compare early mode growth, morphology and refinement; no pointwise exact solution presumed.
4. Circular Sedov/Blast, converging2D Noh, cylindrical obstacle flow using original wall geometry and physical configuration. Measure shock radius/symmetry, central/wall error and wall work/impulse.

Every stage checks convexity, positive volume, exact shared boundaries, fixed per-cell material mass/no matter transfer, total mass and momentum-energy accounting including boundaries/body forces, and physical accuracy at the registered target. Failed convexity must not be hidden by clipping, steering, merging, remapping or mass exchange. Native C/C++ kernels for computation; Python only orchestrates/assesses/plots.

## Implementation boundary

In1D, one contact speed per interface can move a shared point and supply pressure work. In2D, independently assigned face normal speeds need not produce one compatible velocity at a shared polygon vertex. The vertex motion and force/work discretization therefore require a compatible multidimensional coupling. Merely averaging velocities and retaining old face work would recreate the diagnosed mismatch. No claim of exact power/Laguerre representation is made for an arbitrary material polygon. Existing geometry and affine-moment helpers should be reused where applicable, not rewritten or discarded.

Before whole-domain runs, verify actual area derivative from vertex motion equals face-integrated normal velocity, force/work balance and1D reduction on small convex meshes. Those controls are prerequisites, not substitutes for evolved2D tests. Reuse archived internal rotation only with a verified geometry/energy connection; no new empirically tuned stress sensors.

## Parallel manuscript task

One user-requested agent owns manuscript integration of the current1D results, methods and limitations. Root owns2D implementation and tests. Manuscript update must not recast1D evidence as2D or Kepler validation.
