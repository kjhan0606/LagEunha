# Log primitive reconstruction: scalar evidence and inherited velocity blocker

The actual-stencil scalar interval was rejected on a measured finite topology
jump; see `COMMON_STENCIL_PRIMITIVE_2026-11-09.md`. The following two isolated
probes explore a variable choice rather than another physical damping term.
No production, Kepler or boundary executable was replaced.

## Coordinate choice and acceleration relation

For positive q=rho or P, reconstruct ell=log(q/q_ref). The reference merely
fixes units: multiplying q by a positive constant adds a constant to ell and
does not change its LS differences. Reuse the existing face-local tangential
interval and ordered-jump projection in ell. Exponentiate only after projection.
No pressure floor is added; underflow to zero or nonfinite output rejects the
state. This is not the historically rejected mass-exchanging entropy-variable
KEP/ALE correction. The exact Riemann and direct force/volume/work remain the
same, with no mass transfer.

`common_log_primitive_probe_20261109_v1` first kept the old density argument in
pressure-offset geometry. It passes the retained static failure and scalar
topology tests, but does NOT preserve the intended common-acceleration model:
using grad(log P)/rho is different from using grad(P)/rho. It is frozen as an
intermediate diagnostic, not advanced into evolution.

`common_log_acceleration_probe_20261109_v1` transforms the same relation:

    a = grad(P)/rho = grad(log P)/(rho/P).

Thus the existing pressure-offset function receives rho/P for its two pressure
reconstruction arguments. Physical cell density, EOS and force computation do
not change. This is the chain rule for the fitted variable, not a new coefficient.
It does not prove exact reproduction of affine physical pressure at finite
resolution; that earlier balance property must be tested, not inherited by name.
Smooth second-order accuracy and strong-shock positivity over time are unproven.

## Completed scalar and 1D checks

The full static Sedov diagnostic state returns status1 rather than vacuum2,
with total energy-rate residual1.42e-14. This is not a Sedov evolution result.
For fixed positive nonsmooth pressures on the saved vanishing-face geometries,
maximum force difference decreases from1.4959e-5 to1.4961e-7 to1.4961e-9 as
|epsilon| decreases from1e-4 to1e-6 to1e-8. The unweighted-extrema scalar
candidate's finite jump is absent in these fixtures.

All nine 1D endpoint tests finish with exactly fixed cell mass and energy
balance near roundoff. Same physical ICs are byte-identical to the controls;
N=200 except Shu-Osher800, with original exterior guard cells and times.

|Case|Density relative L1|
|---|---:|
|sod|0.003549043283|
|blast|0.03470604581|
|noh|0.004314038745|
|shu|0.003210339493|
|lax|0.003164836114|
|doublefan|0.02300350383|
|collision|0.003771950033|
|contact_open|2.589533526e-13|
|moving_contact_open|1.092093083e-13|

Blast improves from the bounded-velocity parent's.0531785 to.0347060; Sod also
improves, while Noh/Shu/Lax density errors rise slightly. This is mixed accuracy,
not a global method selection. Eighteen actual initial/final states embedded as
weighted rectangular strips agree with the2D interior RHS to
5.192e-14; transverse force is zero. This verifies a dimensional reduction,
not moving2D geometry, boundary conditions, MPI or3D evolution.

## Independent inherited velocity defect

A separate fixed-state test uses rho=P=1 and per-cell velocity components
uniformly sampled in[-.1,.1], held fixed by ID across the SAME vanishing-edge
geometries. Eight seeds are saved. The bounded-velocity parent and the log
candidate give identical results in this test. Their maximum force difference:

|abs(epsilon)|Maximum force difference|
|---|---:|
|0.0001|0.0005733985507|
|1e-06|0.0005731769346|
|1e-08|0.0005731747181|

The limiting difference is8.1024% of the unperturbed maximum cell-force component.
It is inherited from velocity stencil extrema, not caused by log thermodynamics.
Thus the old smooth topology fixture was insufficient for velocity as well as
scalars. A zero-length neighbor must lose its influence on other finite faces.

**Decision: do not adopt or start native2D evolution for this combination.**
The scalar result is useful retained evidence, but the common operator is still
discontinuous. Next examine a continuous face-local velocity reconstruction that
bounds its common mode without resurrecting the rejected unbounded-common-speed
projection. No per-test switch, changed tolerance, new floor, or3D run follows.
The C/CUDA parity certificate of the earlier scalar candidate is NOT assigned
to this new source; it has not yet received a new CUDA qualification.

Evidence is retained under the above directories in plan/build/report,
one_d/candidate_accuracy, dimensional_reduction/report and
continuity_velocity/report JSON plus their source and exact fixtures.
The astronomy draft precedes the second technical paper. Full2D/3D consistency,
Kepler compatibility, and <=2x GIZMO at matched accuracy remain open.
