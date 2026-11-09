# Bounded native periodic weight KDK — 2026-11-06

The existing native C/MPI/CUDA code now has a default-off short smooth
qualification entry: GFS_NATIVE_WEIGHT_KDK_VALIDATE=1. It is restricted to
one MPI rank, GPU variational pressure, periodic boundaries, explicit active
weight state, no viscosity/steering/external force, and a restart window <=.02.
It neither replaces the qualified W0 Kepler branch nor enables shock evolution.

## Method

Half-kick both physical velocity and weight momentum; drift positions and W;
rebuild power geometry; update U by the existing fixed-mass adiabatic map;
recompute force from the refreshed U/V; second half-kick both momenta.
Total energy includes pi^2/(2*inertia).
The initial native .25/omega timestep bound is retained, not loosened.
Later bound violation/geometry failure causes full-payload restore and halving,
up to 12 attempts. No mass or energy repair is applied.

The native periodic state refresh is mandatory after the volume/adiabat update.
An initial engineering execution without this mandatory refresh had larger
energy error; it is preserved as results_v1_pre_refresh.json, not used for
the convergence claim. The corrected v2 kernel refreshes both owned/padding
thermodynamic fields before force evaluation.

## Executed tests

Artifacts: Research/MaterialVolume_20261105/native_weight_kdk_20261106_v1.
Syntax A10; one process, two OpenMP threads; Gresho 16x16 and KH 16x32.
Both advance from t=0 to .002 and then independent timestep-convergence runs
advance from t=0 to .02. This is integrator qualification, not a resolved KH
growth or Gresho accuracy test.

At t=.02, dt caps 1e-4 and 5e-5:
- Gresho relative augmented energy error: 1.527e-11 and 3.818e-12.
- KH relative augmented energy error: 1.274e-13 and 3.194e-14.
- All individual masses exactly unchanged in checkpoint comparison.
- Volume sums equal 1 and 2 to roundoff.
- Maximum per-cell adiabatic entropy change <1.5e-14.
- W and pi become nonzero; this is actual weighted time evolution.
Energy error decreases by approximately four when dt halves.

With the native timestep and no artificial dt cap, both reach .02. Gresho
takes 39 accepted steps; KH 51. Most steps retry once because the later bound
is slightly smaller than the initial bound. These costly rejected trials
must be addressed before claiming any performance improvement. No GIZMO
runtime comparison has been made.

## Limits and next work

Current alias mapping requires ordered dense IDs and uses a scalar box length;
it must be generalized/validated for arbitrary reorderings and rectangular
domains before broader use. The tested KH rectangle has Ly=2Lx.
This is smooth variational pressure with a weight-inertia degree of freedom,
NOT the 1D exact material Riemann/PLM solver. A common shock/entropy formulation
and the physical role of auxiliary inertia remain unresolved.
No KH growth endpoint, RT, 3D, shock pass, or <=2x GIZMO claim is justified.
Preserve all existing 1D and Kepler results.

## Rectangular/permutation update

Artifacts: Research/MaterialVolume_20261105/native_weight_alias_20261106_v1.
Reused existing gpu_unified_owner_aliases ID hashing, independently of the
unified force model. Native periodic refresh now passes Lx and Ly separately
and uploads the certified owner index map. Device state equality checks and
immutable geometry-W checks remain active. KDK saves/restores original array
order but matches state by sorted ID lookup rather than assuming ID=index.

Four short t=.002 evolution checks pass: original KH layout in both forward
and reverse storage order, and y-compressed Ly=1.3 rectangular fixtures in both
orders. The compressed fixtures are boundary/ordering tests, not standard KH
physics results. Final physical-state comparison is in permutation_audit.json.
The prior dense-ID/scalar-box limitation is resolved for these one-rank tests.
Multi-rank aliases, arbitrary domain origins, and full physical endpoints
remain outside the demonstrated scope.

Forward/reverse physical-state differences are <=8.89e-16. Against the prior
same-order build, compared vx/vy/W/pi/U64/V64/mass fields are identical; the
whole checkpoint bytes are not identical, so no full-payload parity is claimed.
