# Native C solve for the linear weight transport correction

Follow-up to KH_VOLUME_TRANSPORT_2026-11-06.md. This is a linear kinematic
and cost screen on its native32x64/64x128 t=1 geometries, not a new PDE run.

## Reuse and equation

Existing source contains CPU-only dynamic face-fitting and volume-Newton
paths in Exam/NativeU64/exam.c. They include historical barrier/target choices
and are not activated. For an isolated linear solve, solve.c includes the
existing Jacobi-PCG implementation from Exam/HydroCore/material_volume_cost.c
with its main renamed. Only its solve/laplacian/reduction routines are called.
The retired psi-times-strain physical construction is NOT used or revived.
Limits remain256 iterations, relative TRUE residual1e-10,5-second cap.

For each neighbor face, L has coefficient area/(2 separation). L is the
actual weight-to-volume-rate graph Laplacian. Given the previous transport
defect R xi, solve

    L q = -R xi,   mean(q)=0.

Only the roundoff mean of the compatible RHS is removed and recorded.
Directed faces are symmetric, and the returned native C solution is checked
independently with a SciPy sparse L. Each face then has its own weight-induced
normal velocity (q_i-q_j)/(2d); equal speeds on all faces are not assumed.

## Accuracy and measured solver cost on Syntax

|cells|threads|iterations|true relative residual|median solve ms|
|---:|---:|---:|---:|---:|
|2048|1|12|2.19e-11|0.465|
|2048|2|12|2.19e-11|0.328|
|8192|1|11|6.57e-11|1.334|
|8192|2|11|6.57e-11|0.856|

Five repeated solves from zero, same fixed matrix/RHS; setup recorded
separately (about0.73ms at8192). This excludes geometry reconstruction,
CPU/GPU transfers, production allocation/interface costs, and extra rebuilds.
It is not an end-to-end GIZMO comparison or a throughput qualification.

To avoid selecting only an easy high-frequency RHS, also test two other
compatible RHS on the same8192-cell graph, normalized to the same norm:

* Long sinusoidal input:90 iterations, true residual9.80e-11, median6.057ms.
* Zero-mean random input: reaches256-iteration cap with true residual4.62e-7,
  fails. No valid solution or successful-solve time is reported for this case.

Thus the cheap11-iteration result is specific to the measured short-wave
defect. General robustness/cost is NOT established. Do not raise the iteration
cap and treat that as satisfying the runtime requirement.

## Admissibility and energy limitations

Rescale the displacement to maximum magnitude1e-4 grid spacing and use a
diagnostic interval1e-4. Linear predicted max |delta W_i-delta W_j|/d^2 is
2.13e-7 at32 and6.12e-8 at64, starting at W=0. This is far below the native
.9 margin but is only a local infinitesimal check: no updated tessellation,
moving separation bound, nonlinear Newton step, or long-time accumulation
has been validated. This is not permission to relax the existing margin.

Work reciprocity P^T Lq=q^T LP is verified to2.20e-15 for nonuniform diagnostic
pressure. It is NOT a total energy pass. In fact P^T Lq is nonzero: changing
volume by Lq while retaining only the old B^T P momentum force would add
unbalanced pressure work. A compatible momentum/work coupling must be derived
before time evolution. The test also enforces cell-volume rates, not an exact
material velocity at every point of every face.

## Decision

Scalar weight differences can cancel this measured linear volume-transport
defect without needing equal face speeds. That part is now demonstrated using
native C. Adoption remains blocked by an unspecified nonlinear physical closure,
its force/work consistency, and robust solve cost; this is a research limitation,
not a request for user approval. No production default or physical parameter
was changed. Continue towards the common1D/2D/3D method;2D is not complete and
the <=2x GIZMO target remains unverified.

Artifacts: Research/MaterialVolume_20261105/kh_weight_transport_solve_20261106_v1/.
Native wrapper, build command/log, inputs, solutions, independent metrics,
and broad_results.json retain the failure rather than hiding it.
