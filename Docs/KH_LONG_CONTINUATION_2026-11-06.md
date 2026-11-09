# KH continuation: frozen completes, dynamic geometry gate stops

Artifact: Research/MaterialVolume_20261105/kh_weight_long_20261106_v1.
16x32, Syntax A10, MPI1/OpenMP2, existing native code. No new force/heat
terms. Resume existing t=.3 states; dynamic first performs a .3->.315
regression and continues that new checkpoint, avoiding its replay.

Changed only the experimental native KH run-duration gate: default .02
remains; explicit GFS_KDK_MAX_WINDOW is now accepted up to2 for this path.
Existing W0 path behavior is unchanged. CFL, later-stage timestep checks,
geometry rejection, mass/entropy checks and 12-attempt rollback remain.
Used max window1.2, endpoints .6,1,1.5. This removes repeated process startup
and endpoint/bin resets; it does not increase the accepted CFL.

Short regression old/new .3->.315 has identical entire particle payload
except the first 8-byte next linked-list pointer, and identical footer.
Raw checkpoint bytes are NOT identical because pointers and checksum differ.
The initial strict comparison correctly stopped; after verifying the actual
struct layout, continue.py checks all non-pointer payload bytes before running.

Frozen completes t=1.5. Segment relative energy changes:
.3->.6 -3.6892e-8; .6->1 -1.6135e-7; 1->1.5 +1.7249e-6.
Cell masses unchanged, segment entropy change <=4.71e-15.
These are not proof of KH growth accuracy or spatial convergence.

Dynamic completes t=1, but fails the next segment. Last accepted state in
log t=1.0731293139821652; retry budget exhausted rc99 after repeatedly
shrinking the step. GPU reports volume=0 velocity_fit=0 geometry=16 halo=0.
Here 16 is a count, NOT a bitmask. The aggregate geometry count includes
weight clamp hits and other geometry failures, so the log cannot yet isolate
which subcondition failed. Do not call this proven cell disappearance.
Last accepted log has Vmin=.0037543517824253805 and Umin=.056511130270904757.

Last PERSISTED accepted dynamic restart remains t=1.0; the restored in-memory
state near1.07313 was not written before abort. No long process remains live.
Raw results.json reports the unchanged input checkpoint after failure;
audited_results.json explicitly nulls its misleading zero change metrics and
records the logged accepted time separately. Do not compare frozen1.5 with
dynamic1.0 as though equal-time results. Preserve all failure logs.

Decision: do not select dynamic weight inertia as common/default hydro.
It is slower and has now failed an actual KH continuation despite conservation
being good at earlier times. This does not reject time-dependent Laguerre
weights in general. Keep frozen W0 as control; isolate geometry failure before
changing weight law, and run higher-resolution KH for physical convergence.
RT, common 1D/2D shock method, 3D and GIZMO timing remain incomplete.

## Terminal status of historical 64x128 PCG-certificate continuation

Reinspection of `material_kh_resolution64_resume_20261106_v1/1.3/log` confirms terminal MPI_Abort99, not a running job. Last stage-start time1.2715245061727007; target1.3 NOT reached. Last trial reports active20, remaining4 atdt4.8828125e-8 after previous changes of active set; the subsequent retry ceiling is2.44140625e-8. The diagnostic failure snapshot explicitly is NOT an accepted restart. Preserve all files; do not restart it or report it as current thermo-only candidate success. Lineage identifies historical binary `material_pcg_certificate_20261106_v1/probe.exe` SHA25645f3ae20c07dcc68a84767ffde71732bc484b8125331fc331a03ff3cc205da63, resumed from accepted t1 checkpoint. This is a separate unresolved constraint/geometry failure from current RT finite-width runs.
