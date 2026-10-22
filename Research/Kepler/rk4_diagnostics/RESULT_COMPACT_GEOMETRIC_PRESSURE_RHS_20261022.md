# Completed-space geometric-pressure RHS and original CUDA/MPI mass solve

2026-10-22. Native implementation continues under Astra's compact transaction
and curved-domain addendum. The preceding turn was progress: geometry/mass
components and their fixtures passed. This turn connects the finite-functional
thermal, volume-adjoint, force and original resident-mass components. No new
Astra call or separate accepted-source guard experiment was made.

## What changed

`material_compact_volume.h` evaluates canonical P3/enriched exterior traces
in the geometry-completed space. It reuses the original endpoint-subtracted
area functional and its analytic G/H derivatives, including both dynamic
geometry modes. Mean-pressure force and external-pressure work use this same
finite boundary functional. The difference from direct physical cofactor
traction is measured; quadrature accuracy is not inferred from compatibility.

Physical exterior is an explicit property of the actual domain registry.
A one-sided LOCAL face may have a cell on another MPI rank. Such a face
receives no external P0 traction. Its local cell contributes its pressure
adjoint; original owned/ghost summation combines the two sides.

`material_compact_thermo.h` uses immutable projected physical internal energy
per parameter-reference volume, e0. For positive corrected volumes j0 and j:

    e_ref = e0 * (j0/j)^(gamma-1)
    P = (gamma-1) * e_ref/j
    de_ref/dt = -P * dj/dt

The generating pressure may differ from P0. No projected heat is reset to
P0*j0/(gamma-1). The original delta-pressure/expm1 decomposition is retained
for cold-pressure precision. Mass density has an explicit reference-measure
tag, avoiding a second Jacobian factor.

`material_compact_thermo_cuda.h/.cu` implement the original three-pass
geometry -> pressure -> residual-pressure/gravity structure using the new
basis and positive rho/e fields. Original MFEvolvedBlock/19-metric ABI and
compensated reduction pattern are reused. The original gfs_material12 stress
operator evaluates the same position/velocity/test space. No viscosity,
steering, independent spin or energy reservoir is introduced.

`material_compact_operator.h` supplies persistent device arrays, exact
owned/ghost state-bit checks, per-cell current delta/delta_dot/mean pressure,
all three CUDA passes, canonical mean-pressure/exterior forces and original
MPI force accumulation. It reports continuous energy/work/torque identities.
Reference fields, generating coefficients and functional remain immutable.
This class has no state-commit, restart or trajectory API. It assumes the
canonical partition and continuous proofs established elsewhere.

## Actual new test evidence

`material_compact_operator_v1_ada_mpi2` is a NEW two-cell/two-rank finite
algebra/API fixture on authorized lageunha. It does not read an accepted
checkpoint, rerun a previous physical interval, project a real source field,
or validate long Kepler evolution. Its curved exterior uses the ORIGINAL
nonpolynomial enriched radial functions with synthetic frozen coefficients.
Unlike the previous polynomial-only fixture, both an enriched curved exterior
and a shared MPI interface are exercised.

The fixture has 18 global owned DOFs and 60 points/pass. It uses nonuniform
positive projected rho/e, arbitrary dynamic geometry velocities, pressure,
central gravity and explicit physical exterior Pext=P0. It compares all19
metrics/forces against a native reference, differentiates the thermal energy
independently, then solves both acceleration components with the preserved
original resident CUDA/MPI mass solver.

Actual result, native0/runner0, wall2.440s:

- CPU/GPU normalized metric difference1.77235e-16.
- Thermal difference/U7.98919e-17; force difference1.65172e-17.
- Continuous RHS energy identity/U4.68380e-13.
- Torque identity2.18414e-18.
- Independent finite thermal derivative error/U8.32292e-14.
- Canonical adjoint/direct-cofactor difference4.96505e-17.
- Maximum absolute delta3.27135e-5 in this deliberately finite functional.
- Original resident mass:18 inner iterations across both components,
  binary80 true residual5.43779e-20.
- After mass solve, energy identity/U4.68476e-13.
- Shared MPI face correctly receives no exterior traction.
- Thermal-seed check confirms e0 is preserved and generating P differs from P0.
- Dynamic coefficient arrays are unchanged by RHS evaluation.

RHS wall1.9383ms is a small-fixture measurement, NOT a production extrapolation
or the amortized t1 cost gate. Positive sampled J and algebraic closure do not
certify whole-map injectivity, continuous thermal positivity or quadrature
accuracy. A delta value in this fixture is not evidence of physical density
accuracy in the real source.

## Builds and preservation

`material_compact_rezone_build_v7` passed native MPI/OpenMP compilation and
new AdaSM89/H200SM90 thermal CUDA compilation. Component links include the
unchanged v5 geometry CUDA objects, earlier sampler objects and preserved
original Ada/H200 GPU objects. Existing GPU APIs and the new batch API resolve.
37 project compiler dependencies, frozen/root source equality and all source,
object, retained-input and fixture manifests pass. This is a relocatable
component link; no complete transaction executable is claimed.

Sessions64716,55668,2973,51339,29407 and5934 are closed0. mpirun used
</dev/null; its native/runner return-code receipts are both present. No new
simulation, Slurm job, background dispatcher or monitoring loop was launched.
Original repository HEAD8e25dfe, old kernels/binaries/source/checkpoints and
shared global context are unchanged.

## Actual physical state and remaining work

The accepted trajectory remains step5 at t=0.00011263377801805954, prefix
`material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r`.
No source import/context proof, actual target partition, conservative remap,
compact commit or new physical time interval was executed in this turn.

Next implement the actual transported generators/current curved power-cell
partition and finite W-linked mass law, canonical event/coverage and positive
overlap functional, rho/e and coupled shared p/L/K projections with independent
field/entropy/error budgets. Then connect global geometry/rank/positivity
guards, generating-data sealing, restart and the complete bounded native
controller. Actual source force-response comparison remains part of the
conditional transition, not something proved by this synthetic fixture.

Before production, also bound evaluation of very SMALL exterior geometry
remainders: the current point implementation subtracts absolute Gamma/chord
and its P3 interpolant. Normalization can amplify cancellation near the real
tiny step5 deformation. A stable relative/compensated trace evaluation and
canonical subinterval endpoint handling may be required; the present large-
amplitude fixture does not establish that bound. This is an implementation
accuracy condition, not permission to flatten the exterior or drop a mode.

Only a passing actual transfer allows the next native RK4 step. Include
actual collar/remap/guard costs in the unchanged amortized t1<=24h gate.
Do not restart old10.8M-space long production or rerun completed preprocessing.
t6/8/10, repeated-remap lifetime, shock/W0/boundary/paper normalization and
the longest GIZMO test remain open. Artificial shear heating is unresolved.
