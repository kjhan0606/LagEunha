# Exact repeated-work reuse in the current common operator

This is a computational candidate, not a new physical model and not an adopted production default. The existing gfs_material_rhs_memo2d.h and gfs_material_weight_memo2d.h are reused; historical qualification is documented in MATERIAL_RT_PORT_AUDIT_2026-11-06.md (pure-RHS memo and combined-weight memo sections). Prior RHS-only timing improved about1%, combined about10%, so neither is presumed to give a new speedup without measurement.

## Scope and build

Only the frozen candidate's gfs_material_trace_probe.inc is reconnected to the existing exact-input caches. Keys contain complete prepared face/entry/cell data, gamma/limiter data for the pure RHS, and edge/cell/RHS/wall data for the weight solve. Hits copy every output of the respective pure function. Padding may cause misses, not false matches. The pressure/volume equations, tolerances, geometry, actual timestep bounds, gravity, work ledgers and endpoint acceptance still execute through the same code.

Both GFS_RHS_MEMO and GFS_WEIGHT_MEMO require the explicit value1; absent/0 stays off. RHS response-export calls bypass memo when an export buffer is required. Stateful GFS_MATERIAL_WEIGHT_RECYCLE=1 bypasses weight memo to preserve numerical-history side effects. Only serialized single-rank use with no device changes inside the process is in scope. General MPI/thread/3D reuse is not claimed.

Build: baseline SHA256 11a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e; candidate 9602d021ffcd89456cdec12f0a671ecb03d7a0d62b0014b8730a135ade9b2037. Exactly one archive member changed, exam_gpu_extract.o; CUDA objects and libcyl remain unchanged. Candidate.patch/build.py/build_audit.json retain the source/build lineage. Canonical source, live Cylinder/KH/RT executables and configurations are untouched.

Existing CPU key-mutation fixtures were compiled against this candidate's frozen headers and passed. They test exact repeat hits and invalidation for state/geometry/gamma and RHS/geometry/velocity changes, comparing actual pure solver outputs. These CPU checks do not qualify the native CUDA evolution.

## Ordered H200 qualification, pending

Job456838 requests one H200,2 CPUs,8GiB,15min and depends on successful job456829. This dependency is necessary: it consumes the current portable H200 reference fields and requires its six architecture gates to pass. Latest scheduler state is PENDING/Dependency, afterok:456829 unfulfilled. Job456829 is pending QOSGrpGRES. Neither has been reported as executed.

The finite job performs20 native runs sequentially: default-off and enabled gates for KH,RT,Sedov,fixed-solid Cylinder and Noh; RT64 and RT128 same-state ABBA windows; and a recycle-enabled bypass comparison. All13 saved fields and accepted time/dt sequences must match bitwise within each comparison, with exact masses, positive U/V, correct endpoint and volume error<=1e-11. Failure stops the batch with evidence retained. RT128 uses the accepted initial preflight checkpoint at.002, not a rewrite/restart of its live full history. Twenty new final checkpoints plus copied inputs are budgeted below1GiB. All paths are /gpfs; no /scratch use.

The timing summary reports whole accepted-step wall, excluding the first warm-up step, and cache hits. More hits alone are insufficient for adoption. Shared-node/short-run results would still not prove whole-trajectory physical accuracy or <=2x GIZMO at matched accuracy. H200 results, failure/retry coverage and meaningful total runtime improvement are required before promotion. The second paper remains deferred until the first astronomy draft.
