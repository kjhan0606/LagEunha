# Actual-stencil primitive projection: rejected before 2D evolution

The preceding goal turn was evidence progress: actual glass IC generation and
matched regular/glass outcomes were independently verified. This turn completes
the pending candidate tests and supplies a new counterexample that changes the
next action. The full common-method goal remains active and incomplete.

## Candidate and scope

Evidence directory: `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_stencil_primitive_probe_20261109_v1`.

Replace the scalar-v3 tangentially shifted scalar interval with each owner's
actual reconstruction-stencil minimum/maximum. Apply the existing bounded-pair
projection to density/pressure as scalar components, as already done to normal
velocity. Unlimited face-measure LS, pressure-offset fit, exact Riemann solver,
force/volume/work exchange and SSPRK2 stay unchanged. No new coefficient, floor,
viscosity, test-specific switch or mass exchange was introduced. The old
unlimited LS is inherited, not a newly removed gradient limiter.

The periodic C/CUDA research sources are isolated. Production, existing Kepler,
wall/RT executables and running states were not replaced. No native 2D evolution
or 3D run was started for this candidate.

## What passed

- The retained Sedov failure state now evaluates with finite RHS and no vacuum
  status. Total-energy-rate residual is 4.15e-15. The fourth output column is
  H=Udot+P Vdot, not total energy; its sum 88.42 is not an energy defect.
- Previously recorded smooth wave/disappearing-face fixtures remain continuous.
- Nine 1D cases reach their original endpoint with exactly fixed cell masses;
  maximum absolute energy balance error is below 4e-15. Accuracy is mixed:
  - sod: density relative L1 0.004007015214.
  - blast: density relative L1 0.05242686441.
  - noh: density relative L1 0.00422554947.
  - shu: density relative L1 0.002971852611.
  - lax: density relative L1 0.003028359178.
  - doublefan: density relative L1 0.02312788621.
  - collision: density relative L1 0.004076504176.
  - contact_open: density relative L1 2.589533526e-13.
  - moving_contact_open: density relative L1 1.071531752e-13.

- Twelve full-RHS C/CUDA fixtures agree within 1.83e-15 scaled error.
- Eighteen initial/final 1D states embedded as weighted rectangular 2D strips
  agree in interior force, volume rate and total-energy rate to
  4.582e-14; transverse force is zero. This excludes outer boundary stencil
  differences and does not certify global moving geometry, MPI, or 3D.

## Falsifying geometry test

`continuity_nonsmooth/run.py` holds rho=1, velocity=0 and positive pressure
exp(uniform(-1,1)) fixed by cell ID for seeds0..7. It reuses exactly the existing
1024-cell periodic shear geometries with epsilon=0 and +/-1e-4,1e-6,1e-8.
Only geometry changes; no evolved solution or physical coefficient is changed.

Maximum force difference from epsilon=0, over both signs/eight seeds:

|abs(epsilon)|Scalar-v3 parent|Actual-stencil candidate|
|---|---:|---:|
|0.0001|1.377926177e-05|0.01214597359|
|1e-06|1.378047201e-07|0.01214994079|
|1e-08|1.378048418e-09|0.01214998046|

For seed0 and epsilon=-1e-8, the candidate force jump is 22.8129% of the
maximum unperturbed cell-force component. It does not shrink with the geometric
perturbation. The parent error decreases by100 when epsilon decreases by100.
Thus this is not merely a large gradient or a roundoff difference.

The actual-stencil extrema ignore face measure. A neighbor with vanishing face
length can still set a finite bound on other, surviving faces; at the topology
change that bound disappears. Smooth-wave fixtures did not activate this defect.
Affine exactness is also only conditional on face values lying in both stencil
ranges; the passed rectangular-strip test is narrower than arbitrary geometry.

**Decision: reject this candidate, preserve all negative evidence, and do not
expand it into native 2D time evolution.** Positivity at one failed state and
CPU/CUDA agreement do not compensate for a discontinuous physical RHS. No
per-test fallback or coefficient adjustment follows. The next design must address
positive scalar reconstruction without an unweighted topological min/max.

Second-order convergence, common boundary treatment, the full2D test matrix,
3D verification (after2D), and <=2x GIZMO at matched accuracy are not achieved.
The first astronomy draft remains prior to the proposed second technical paper.
