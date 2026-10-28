# Native preimpact fixture and local wall-RHS extension

The existing Q2 polygon-wall failure has been reproduced with its **unchanged native CUDA+MPI+OpenMP executable**, starting from the saved t=47.5003165356 state. A 438.419-second Ada replay captured the full accepted predecessor at **t=48.37534965270793**, step 11580. The next trial hits the original wall guard for the same generator **55283**, facet **543**, with the identical recorded crossing coordinates. The accepted predecessor SHA-256 is `c4b11f1996349f5ce8ab664569fa4415cd91b5f836fe9579c98f7dcd51aff758`. No t0 or t30 restart was required. No collision repair has yet been applied.

Fixture: `/gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/control_capture/khout.011580.dat`.

## Why a trial RHS extension is possible

For fixed domain Omega, the cell is `Omega intersection Voronoi_i(q)`. A generator slightly outside Omega can still own a positive clipped cell. An independent 16-cell calculation checks just-inside, on-plane and just-outside configurations. The force matches minus the energy gradient to about 8.3e-10 at the smallest finite-difference interval; static-wall power closes within 2.8e-17. This is a local energy/geometry result, not an assertion that arbitrary external positions remain valid.

The original native GPU code rejects contact/outside sites before its first RHS because two clipping loops explicitly require positive generator-to-wall plane distance. Original probes reproduce that rejection. A separate candidate changes those checks only under a **K1-only diagnostic mode**, with fixed maximum excursion **1e-6**. It serializes first-RHS state and exits before any predictor or accepted evolution. Band overflow, missing/empty cell, nonpositive state, incomplete support and capacity guards remain failures. The candidate binary SHA is `3a9e98eaed65a903534a001f6c400c4ad815b7d873577ad52cd658ab1b580e45`.

Astra revised its recommendation to a restricted trial-only extension followed by actual RK trial event bracketing. This avoids treating a predictor crossing as a physical collision. It does not yet approve or validate an implemented event controller; rollback and grazing/event accuracy tests remain outstanding.

## Actual native probe results

Five finite first-RHS configurations pass: original predecessor, 1e-7 inside, contact, 1e-7 outside and the configured band maximum minus a 1e-12 representability margin. All 55351 IDs, coordinates and masses are retained and physical fields remain finite with positive U and volume. Independent bounded polygons for the crossing site and its neighbors match native volume within **5.47e-8 relative**. Independent pressure acceleration differs by at most **4.03e-19**, and thermal work by **3.45e-15**.

Selected native face records cover four owners and **74 support keys** in each configuration, matching all independent keys. Maximum absolute length/vector/first-moment error is **1.97e-13**; reciprocal-face error is at most **1.06e-15**. Native face-work formulas close within **7.95e-22**. Capturing these records leaves the complete first-RHS state bitwise unchanged apart from the excluded next pointer. These independent geometric tolerances are not a relaxation of the separate failed halo20 bitwise gate.

Two negative fixtures also pass their rejection tests: excursion above 1e-6, and a deliberately empty clipped cell created by two nearby external diagnostic sites. This checks rejection before RHS publication, not yet transactional rollback.

With both new settings absent, the new executable reproduces the old tiny accepted replay's **entire stored record bodies and all 52 footer bytes exactly**, excluding only the linked-list next pointer. This checks ordinary-path preservation at this late fixture; it is not a long-time event result.

CPU independent audits **1236226** and **1236471** each completed 0:0 in two seconds. An initial one-second postprocessor failed because the established snapshot reader does not expose RHS fields; its replacement uses the already captured native stage TSV. GPU probes and physical trajectories were not repeated for that reader correction.

## Output correction and preservation

The selected-face debug switch also enabled an existing all-candidate structural dump by default. That auxiliary capture was stopped after this was observed; only its verified runner/descendants were signaled. Completed baseline/inside face captures were reused. The three remaining captures used the already available `GFS_GPU_TESS_DUMP=0` option. Exactly the three newly generated, unnecessary all-candidate text logs were removed (5,412,487,944 bytes); selected face records and all native states/inputs were retained. No physical simulation or unrelated process was stopped. The cleanup receipt lists every path and size.

## Next step and remaining goal

`KEPLER_SLIP_WALL_NATIVE_2026-10-29/EVENT_CONTROLLER_DESIGN.md` records how to add accepted-time RK4 event search, full-payload/cache/ledger restoration, no predictor reflection, first-outward-crossing detection including out-and-back paths, bounded retries and separate wall impulse/torque accounting. The frozen diagnostic worktree and binary must be preserved; develop event integration separately.

The slip impulse header was also made C-compatible and revalidated in C, C++ and 4096 actual CUDA cases. This is an impulse primitive, not a complete wall evolution.

The failed periodic Q1 t=3 checkpoint remains preserved and its radial gate remains failed. No t4 production was submitted or acceleration promoted. A completed Q2 event test will validate collision handling only; the next physical campaign must preserve the original Q1 disk in a consistent boundary IC and pass t3/6.5/10, boundary sensitivity and the matched GIZMO t120/t600 plus interpreted longest-orbit target. Long-time stability remains unachieved.
