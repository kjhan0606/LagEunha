# Common search-capacity correction — 2026-11-09

The preceding implementation goal turn made evidence progress through the reconstruction ablation and qualified RT recovery. The intervening glass discussion verified existing results. This turn changes authoritative source and completes actual regression tests; the full common-method goal remains active.

## Adopted change

The previously diagnosed RT failure was a conservative search-radius overflow, not an invalid physical cell. The existing native geometry transaction now consumes the GPU's measured pure-halo shortfall and extends its computational grid and periodic/reflected padding together. The same width is used by the force halo refresh. The same fluid stage is reevaluated before volume writeback; mass, velocity, weight, timestep and physical boundary laws are not modified by this retry.

The existing integer-bin helper is reused. Uniform growth is required because the existing image callbacks accept one scalar width. Context/domain/bin/mode/time-rewind changes invalidate cached capacity. Capacity persists between calls of the same material driver but is inactive outside it. Effective bounds are also carried into force fallback reconstruction and endpoint evaluation. Cylinder keeps its existing callback; finite-exterior Noh does not acquire fictitious extra samples. Mixed invalid geometry, no growth, overflow, unsupported rank count and demand beyond the existing single image layer remain failures.

Authoritative files:
- `Exam/NativeU64/exam.c`
- `Exam/NativeU64/gfs_material_ssprk2_native.inc`

Only archived `exam.o` changes; all CUDA objects and the pressure/force/work, reconstruction and volume inverse laws remain identical to the frozen portable control. The new executable is under `common_search_capacity_20261109_v2/native/native.exe`, SHA256 `cdbd4824693259bde44ee447fb34e35d6d39b57459e76e057fe1bfcad0d81d53`. The ongoing RT process still uses the original frozen executable. No process was stopped or restarted, and no executable used by a running process was overwritten.

## Actual regression evidence

`qualification/results.json` contains nine successful runs:

1. Actual failed RT step, t=3.095874084472433 to3.0959365844724327. One measured deficit3.759e-7 causes halo0.015625 to0.0234375. All13 saved state fields are bitwise equal to the previously verified fixed-halo3 result. Volume mismatch2.9374e-12 passes the unchanged1e-11 threshold.
2. Existing KH, Sedov, RT, Noh endpoint, Noh stage, and cylinder short regression cases. All six reproduce their original endpoint fields bitwise.
3. A nine-accepted-step RT pair, fixed-halo3 original versus automatic-growth candidate. Final states are bitwise equal, with only one capacity extension across the entire automatic run.

`context_audit/report.json` tests actual candidate routines for scope activation/deactivation, next-step retention, context/mode/domain/bin/time reset, invalid deficit and overflow rejection without mutation, single-layer limitation, and cylinder delegation. The preexisting generic integer-bin helper is exercised in dimensions1/2/3; this is not hydrodynamic evolution in those dimensions.

`periodic_growth` uses initial halo1 on KH/Sedov. Both reproduce the baseline, but neither naturally activates growth. That gap is explicitly covered by `periodic_forced_rebuild`: an isolated TEST executable discards its first otherwise valid geometry before writeback and requests one extra bin. Both complete grid/padding/force-refresh rebuilds and preserve all endpoint fields bitwise. This fault-injection source/executable is not adopted. Its artificial shortfall is distinguished from the actual RT failure.

The initial v1 cylinder test omitted the required `.exterior` checkpoint sidecar and stopped during input loading, before fluid evolution. That failed harness record is preserved. The corrected v2 copies both `.stream` and `.exterior`; the cylinder result then matches the archived result. It was not repaired by changing hydrodynamics.

All tests are finite, sequential Ada runs using the existing C/OpenMP/CUDA program. Each case records command, effective parameters, environment, child PID, input lineage and terminal result. The output bound was32 checkpoints/335MB per recorded batch with approximately149TB free; actual batches wrote fewer. Existing simulation outputs remain unchanged. Adoption hashes and source backups are in `adoption.json` and `canonical_before`.

## Live calculation and limits

At14:24:24 KST, supervisor958973 and native959091 were verified alive via/proc. The original RT128x256 trajectory had accepted t=3.2754365844720144 with dt6.25e-5 and zero retries; current segment target3.3, batch target4.0. This is a timestamped observation, not a continuing live-state assertion.

This correction removes an avoidable computational-capacity failure without altering the physical operator. It does not qualify multi-rank MPI, complete2D morphology/accuracy,3D hydrodynamics, or <=2x GIZMO at matched accuracy. 3D evolution remains deferred until2D completion. The rejected log-pressure reconstruction remains unadopted; see `COMMON_ADMISSIBILITY_ANALYSIS_2026-11-09.md` for the new exact-solver counterexamples before another reconstruction change. The first astronomy draft remains prior to the second technical paper; these audit artifacts preserve implementation evidence for the latter.
