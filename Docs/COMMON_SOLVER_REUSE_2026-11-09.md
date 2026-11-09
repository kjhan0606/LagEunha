# Reuse of a numerical volume solve

The private native candidate reduces repeated factorization of the existing
volume Jacobian. It changes no physical coefficient, state variable, pressure
force, boundary law or thermal update. The ordinary20-step comparison and
actual inflow/outflow comparison retain all13 extracted endpoint fields
exactly. Near-contact comparisons differ by at most4.53e-9 in density and
4.50e-9 in pressure. This is numerical solver qualification, not completion
of Cylinder accuracy or a common1D/2D/3D method.

## Justification, necessity and conditions

The same full anchored equation J deltaW = Vtarget - Vactual is solved, with
face coefficient area/(2*distance) and the current raw face-difference
residual. Repeating symbolic/numerical factorization dominates part of the
existing cost. Reusing a previous LU as a right preconditioner removes some
repeated algebra without adding a physical model or fitting parameter.

Reuse requires identical ordered material IDs, exterior-anchor rows, and
positive assembled CSR structure. All are measured from the current actual
face/particle data. A changed basis or graph uses a fresh current factorization.
The previous LU is never treated as the current Jacobian. At most8 GMRES
iterations are attempted; this is a computational work cap, not a physical
coefficient or relaxed acceptance tolerance. Nonconvergence falls back to
current factorization; an old-factor error clears the cache before rebuilding.
Every successful solve still meets the existing current raw-J residual<=1e-10.
Existing nonlinear volume, stored-weight precision and physical checks remain.

The first source-only prototype gated only IDs/anchor rows. Review found a
counterexample: changing internal connectivity can isolate an unanchored
component while a compatible RHS still has a small residual. That prototype
was not used for native evolution. The qualified version also requires the
same positive CSR pattern; its singular-disconnection unit case is rejected
by the fresh solver's pivot test. This is not a newly observed hydro failure.

A successful Krylov solve has no current PARDISO factorization and hence no
current-factor pivot diagnostic. The current residual and nonlinear geometry
checks are the acceptance evidence. Do not claim identical numeric pivot
behaviour for arbitrary changes in coefficients or matched field bits in all
flows. Solver mode1 logs mean successful preconditioning; mode2 means an
attempt followed by fresh factorization; mode0 is direct. Certificate time
inside GMRES is included in solve, not the separate certificate timing field.

## Executed checks

- Sixteen synthetic/helper cases pass, including same-N ID replacement,
  changed ID order/anchors, changed dimension/pattern, invalid input rejection,
  factor/basis lifetime, eight-iteration fallback, injected old-factor failure
  and singular unanchored-graph rejection.
- With reuse disabled, the rebuilt native binary reproduces the previous
  near-contact17-step result: all13 extracted fields and accepted/retry times
  are exact. Only libexam.a:exam_gpu_extract.o changed; all other archive
  members, including CUDA physics/geometry objects, were checked unchanged.
- Enabled endpoint and stages modes both finish the same interval with exact
  retained masses and IDs, positive U/V, raw-J residual<=8.05e-11, maximum
  volume discrepancy<=3.48e-9 and energy-ledger closure<8e-15. Their accepted
  times are slightly different, and their fields are not bit-exact.
- Ordinary endpoint evolution from3.00025 to3.01025 takes20 accepted steps.
  Enabled/disabled final13 fields are exact. The enabled solver reuses20 of40
  linear calls; all current residuals pass and volume discrepancy<=1.01e-12.
- The actual inflow/outflow interval3.01025->3.0135 takes7 steps. At3.0125,
  160 cells enter and66 leave (N126696->126790); at3.013,94 more leave and N
  returns to126696. Three basis preparations are logged, and enabled/disabled
  final13 fields and event records are exact. This verifies actual ownership
  cache invalidation; constant-N ID replacement is independently covered by
  the unit case, not isolated by this native interval.

Measured accumulated linear-solver work:

| Same-source interval | Reuse off | Reuse on | Reduction |
|---|---:|---:|---:|
| Near-contact17 steps, endpoint |10.418s|6.267s|39.8%|
| Ordinary20 steps |16.463s|13.421s|18.5%|
| Ownership event7 steps |5.515s|4.682s|15.1%|

Ordinary whole-run wall time was53.447s versus49.695s. These jobs shared
hardware with the existing main continuation and other work. The timing is
useful implementation evidence, not an uncontended end-to-end benchmark or
proof of the <=2x GIZMO objective. Matrix analysis still occurs when the
positive graph changes; do not silently replace this gate to improve a score.

## Existing GIZMO evidence and remaining requirement

There is an actual earlier2D KH64x128 comparison through t=1.35, not only a
Kepler estimate: native154.551s, GIZMO107.665s, ratio1.4355. The input fields
match, but GIZMO initializes density with its kernel, native uses five restart
segments, hardware is shared, and equal accuracy was not established. Endpoint
vy-RMS differs60.5%; even the two128-resolution results differ14.6%. Neither
is thereby proved more accurate. An earlier KH128x256 t=.5 ratio1.9641 has
similar scope restrictions. Existing KH/RT binaries and comparison machinery
can be reused after the current common method is qualified.

Read-only squeue/sacct queries during this turn found old H200 benchmark job
454560 still PENDING(QOSGrpGRES), without a start/end time. Its saved script
uses an earlier native binary. The job was not changed and is not a completed
benchmark of the present Cylinder or LU candidate.

The remaining goal is matching physical accuracy, IC/BC and final time with
the selected current method, then measuring full initialization, geometry,
weight solves, hydro, retries/restarts and outputs under controlled resources.
Neither a helper speedup nor a cell-count-only ratio satisfies that goal.

## Source and scope

Evidence is under Research/MaterialVolume_20261105/
material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/
wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/
contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/
volume_coordinate_v1/precision_stop_v1/lu_native_v2:

- native/source/Exam/gfs_material_volume_inverse.h
- native/source/Exam/gfs_volume_lu_precondition.h
- native/source/Exam/NativeU64/gfs_material_volume_inverse_native.inc
- native/completed.json, basis_cases_v1/completed.json
- evolution_completed.json, representative_t301025_v1/completed.json
- ownership_event_v1/completed.json, assessment_v1/summary.json

The private binary SHA256 is
3d91ed870afae2f83116920e610b1b21871f982ae548c93e444a9e73a7df5d63.
No canonical default or live continuation binary was changed. This candidate
is suitable for the next saved-checkpoint performance comparison, with the
same physical-method gates. The coefficient algebra is dimension-independent,
but the native inverse adapter remains single-rank2D Cylinder. Periodic gauge,
general2D, MPI and3D evolution require their own existing planned integration
and qualification; no3D evolution is started before2D qualification.

The first astronomy-paper draft retains priority. All evidence here supports
its validation and prepares the later fluid-mechanics paper; it does not
begin a second manuscript or a new independent test campaign.
