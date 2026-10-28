# Next native RK4 event integration

Status: implementation notes for the next work bundle. The current built candidate is K1 diagnostic-only. No RK4 collision controller, rollback event path or production boundary switch has yet been implemented.

Astra revised its recommendation after the independent 16-cell variational derivative check and native local-extension probes: use a strictly limited trial-only RHS extension, rather than reflecting Euler predictors or guessing non-grazing behavior from sampled acceleration. All cells must remain positive, complete and within a fixed tested excursion band. A site outside the fixed polygon can still have a positive intersection cell; contact itself is not a singularity of the clipped Voronoi energy. The old CUDA implementation explicitly assumed positive wall distance in two clip loops; the diagnostic patch changes those checks only in a K1-only mode with max band 1e-6 and exits before any predictor/accepted evolution.

## Existing source to retain

The frozen diagnostic worktree is `/gpfs/kjhan/LagForce/Kepler_SlipWall_EventsBuild_20261029/project`, native binary SHA `3a9e98eaed65a903534a001f6c400c4ad815b7d873577ad52cd658ab1b580e45`. Preserve it and its linked probe inputs. Build the event implementation in a new worktree using the reconstruction/patch receipt, rather than overwriting the diagnostic binary.

`Exam/exam.c` contains the RK4 core and all stage/final bookkeeping. `gfs_ssp_restore` in `gfs_ssprk3.inc` restores complete native payload/counts/dtold/floor and pair counters, invalidates GPU state and discards padding. The existing RK4 transaction cannot just be enabled: it requires periodic postStage and C_geom=.2, whereas this bounded campaign uses C_geom=.1 and must retain its physical contract.

The current `fixed_polygon_postStage_blend` aborts on predictor crossings. The GPU `gpu_geometry_trial_begin/end/failed` hooks can report trial geometry errors instead of fatal exits. Ordinary flag-off behavior must stay bitwise identical.

## Proposed integration shape

Keep a dedicated default-off bounded-W0 one-rank event state, with a complete accepted-start payload backup and immutable ID mapping. Reject unsupported viscosity, steering, dynamic weights, thermal64/entropy transactions and alternate integrators. The initial boundary is the same 1024-gon; use facet normals and record wall torque. Circle reflection is not valid with polygon clipping.

At each tentative RK4 attempt:

1. Save/restore native payload, local/global counts, dtold, stage floor/pair and energy-audit scratch ledgers. The backup is an accepted state at the driver time. Every rejected trial must restore them exactly and invalidate all device geometry/gradient/padding caches.
2. Permit predictor coordinates only inside the tested fixed outward band. Do not clamp/project/reflect them. A band exceedance or incomplete/nonpositive geometry rejects the trial and shrinks time, rather than increasing the band.
3. Intercept the trial after final combination and final bounded geometry, **before** permanent `[RK4E]` cumulative ledgers, alpha update or accepted-step counters. This permits root search without committing rejected trajectories. The accepted endpoint itself must meet the much tighter contact/domain tolerance.
4. Search for the **first** outward crossing across all sites and facets, including out-and-back paths. Endpoint signs alone are insufficient. A cubic continuous RK4 interpolant can propose subdivisions/root candidates, but actual RK4 reevaluations must establish the event bracket. Its third-order dense accuracy is not proof of fourth-order event accuracy. Unresolved grazing or corners fail closed.
5. Record the actual returned RK4 duration at every probe; a CFL-shortened trial is evaluated at that actual duration, never at its requested cap. Reevaluate competing sites at every trial. Stop on missing time progress or bounded retry limits.
6. At a verified accepted contact, apply the tested plane impulse once, with no position or thermal-energy repair. The next step must rebuild RHS/CFL from reflected velocity. Return the actual accepted physical prefix to the existing KH driver so driver time/checkpoint context stays authoritative; avoid an internal macro interval whose partially accepted time cannot be saved correctly.
7. Keep continuous pressure traction and collision impulse/torque separate. Static elastic wall work is zero. Accumulate and test fluid+wall momentum/angular momentum, total energy and roundoff; never book a wall collision as heat injection.

If root search cannot finish, preserve the last accepted start with its exact physical time; never serialize an outside trial as a restart. If an ambiguous path is advanced through a safe inside prefix, that prefix must go through normal accepted bookkeeping and return to the driver before another event search.

## Required executable tests

- Oblique constant-velocity event: collision time, remaining flight, kinetic energy and impulse.
- Accelerated event: h, h/2, h/4 and spatial/time tolerance convergence.
- Non-colliding grazing/turn-back where a predictor exceeds the wall: no spurious impulse.
- Repeated impacts and injected rejection: no duplicate impulse or ledger contamination; restored bytes/caches exact; progress limits fail closed.
- Frozen native Q2 preimpact fixture: `/gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/control_capture/khout.011580.dat`, SHA `c4b11f1996349f5ce8ab664569fa4415cd91b5f836fe9579c98f7dcd51aff758`, t=48.37534965270793. The unchanged old executable reproduces ID55283/side543 crossing immediately afterward. This avoids replay from t0 or t30.
- New Q1 boundary IC preserving the original disk, then unchanged t=3/6.5/10 gates; boundary-radius/shape convergence and GIZMO t120/t600 plus the interpreted longest-orbit target remain required. A Q2 collision repair does not substitute for these.
