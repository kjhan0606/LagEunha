# t=8 Kepler diagnosis and remedies — bundled Astra review, 2026-10-20

The user asks: first develop our own proposed remedies for the t=8 findings, then ask Astra. This request is that review. Read-only review only; no launches, changes, scans, or extra agents.

## Scope and authoritative state

Current goal: long-time cold Kepler stability beyond the former t=6–10 failure, ultimately the longest GIZMO-paper evolution; preserve the user's own CUDA+MPI implementation, improve efficiency. No generator steering; Re=infinity (no physical shear viscosity). Variable Laguerre weights are allowed. User-approved reduction: W=0 => Voronoi and fixed per-cell mass. Genuine intrinsic spin is allowed, but must include torque/energy consistently. Fable cancelled. One bundled Astra request, not frequent calls.

Canonical source: `/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src` (git origin kjhan0606/LagEunha). This is the active native CUDA+MPI branch. Your earlier full-moment/native-face short-time review was a different experimental path; do not assume those terms are in this binary. Current baseline has no evolving intrinsic-spin degree of freedom. `GFS_LAGUERRE_ROTATION` is only a face-velocity correction.

Active production: Slurm 429482, H200, 4 MPI ranks sharing 1 GPU, old binary preserved; t≈8.5 and continuing toward 10. Run root `/gpfs/kjhan/LagForce/Kepler_AdaptiveStencil_CleanLong_H200_Halo8_InitBounds_20261020/production`. Do not interfere with it. Baseline: 256² across 24²; center12,12; Rin2,Rout8; GM1; Plummer eps0.047; P=1e-6; gamma≈5/3; rho≈1/R smoothed annulus plus1e-3 floor (not Hopkins flat profile); fixed masses, W=0, no steering or physical viscosity; RK4; orbital eta0.06. Inner orbital period17.77889; t8 is0.45 inner-edge orbits, while central softened period0.064. Grid spacing0.09375, eps/h=0.5013.

GPU adaptive stencil and init halo-bounds bugs were already fixed in this baseline. A missing-neighbor explanation cannot be assumed afresh.

## Measured t=8 findings

Audit JSON in this directory: `baseline_t8_true_area.json`. Snapshot: production/khout.013378.dat, t8.00013686, full stride576,65536 cells. Global independent scipy Voronoi audit, restricted to64240 cells with a proved periodic-completeness bound. Stored areas agree to~6e-8 (float precision), including central gas. Stored rho agrees with m/V to~6e-8. Thus t8 density errors are not a plotting or incomplete-cell-area artifact at that instant.

- R<0.5:60 cells, L1 relative azimuthal-speed error0.54479, rms(vr/vc)0.30922, median P0.0018396, peak P0.4150934 (initial P1e-6).
- 0.5<R<1:92 cells, L1v0.28181, rms vr/vc0.21926.
- 1<R<2:1284 cells, L1v0.05076, rms vr/vc0.09335.
- 2<R<4:4284 cells, L1v1.8621e-5, rms vr/vc8.5778e-5; 4<R<8 L1v3.2159e-7.
- Total internal energy0.00531127 versus initial0.000864. Inside R2 has U≈0.00436414 (~82% of U); its initial constant-P U≈1.885e-5. Global relative total energy error≈9.3e-5; cumulative negative-U floor injection≈7.24e-4. Small global energy error hides large local conversion of orbital energy into heat.
- Peak pressure at R0.027957, ID33156, V7.64473e-5. GPU stage volume floor=1e-3*(4h)^2=1.40625e-4, so this cell's actual volume is0.544 of the GPU floor. At the t8 snapshot only one R<0.5 cell lies below that floor. No such cells in0.5<R<8 at that instant. Thus the floor creates an inconsistent stage/final geometry in some cells, but is not proven to explain all core heating.

## Newly isolated implementation defect and work in progress

`voronoi_tessellate_kernel` reads neighbor pden_in and writes out_den; formerly both arguments were ctx->d_parts.den. This is a cross-thread read/write race and violates restrict aliasing. It can perturb density gradients/limiters and face fluxes.

Same EXACT t6 snapshot restarted for one RK4 step twice, 1rank native code: all392739 first-force face keys, coordinates, normals and udotds match bitwise, but pstar differs up to2.4884e-7, max dv~7.8244e-7; see `compare_late_repeat.json`. This proves nondeterministic thermodynamic/flux evaluation in the old binary, not its role as the whole instability cause.

Candidate-array optimization preserves all first-force geometry bitwise too (`compare_late_stream.json`) but changes thread scheduling and hence these races. Its compiled thread stack dropped20384→1952 bytes; shared-allocation wall-time speedup was NOT demonstrated (241s old vs259s new, median tess20.38 vs20.52ms), so do not claim speedup. CUDA+MPI and physics model are retained.

Immutable density-input buffer fix is built as `density_snapshot/eunha2_stable_input.exe` and being checked with two identical t6 one-step restarts (`fixed_late_a`, `fixed_late_b`). The full production binary is unchanged. Source diffs and frozen copies are in `density_snapshot/`. Correction freezes exactly the former uploaded density time level; it does not yet improve stage centering or MPI halo state. GPU native end-step volume update still runs on CPU and dominates 1rank cost (~0.8s versus ~4*0.02s GPU tess). Porting that update is a separate candidate, requiring raw-volume/floor equivalence checks.

## Driver's own assessment and proposed route, for critique

1. Finish and verify deterministic GPU stage inputs first; test cross-rank stage/halo consistency. Geometry is already correct at t8; do not keep repeating geometry fixes or adding viscosity. Treat this as correctness, not proof of Kepler stability.
2. Separate underresolved central gravity from hydro closure: eps is only half an initial cell spacing, central orbit frequency~98 while the disk edge frequency~0.353. Use gravity-only source-error/centroid/second-moment controls and a resolved-core IC sensitivity study without changing the annulus, mass, P, and force law. Local initial core refinement (not runtime steering) could be cheaper than whole-box refinement, but must initialize m=rho*actual V and maintain CUDA search completeness. Changing eps should be an explicitly different-potential sensitivity test, not a hidden cure or proof of the same benchmark. The GIZMO paper does not give a numerical eps.
3. Instrument the actual CUDA work/volume contract: per-cell stage geometric Vdot, work-face volume rate, entropy/P response, radius, raw-vs-floored V, and gravity/hydro energy transfer. In uniform-P smooth flow, mismatch R=Vdot_geom−Vdot_work implies dlnK/dt=(gamma−1)R/V for fixed mass. Current GPU work uses Riemann normal speed on ordinary faces, while generator geometry determines thermodynamic V. Need measure whether this is the dominant growth driver; `GFS_VOL_AUDIT` is CPU-only today, cannot just enable it on GPU.
4. If geometry/work mismatch survives correctness and resolution controls, repair the discrete operator/energy contract before extending to long time. Merely enabling geometric face work at fixed mass is not automatically a consistent ALE scheme; mass flux and volume law must also match. Prior Oct8–9 exact-material-volume Laguerre path hit weight-realizability/validity limits and should not be repeated without a new formulation.
5. Rotation is a possible higher-moment closure, not a damping knob. For cell-centered intrinsic angular momentum S, pair forces applied at shared face centroid give torque (cf−xi)×F and opposite partner torque; orbital+intrinsic L can then conserve. Spin energy S²/(2I) needs the I-dot term and matched face work. At Re∞ no viscous spin relaxation may be invented. In a local pure Kepler field, antisymmetric gradient spin Ω/4 accompanies symmetric strain3Ω/4, so spin-only closure omits the dominant subcell shear (9× its energy for isotropic moments). Critique this inference, and whether your existing full-moment work provides a useful conservative CUDA-portable route rather than a new Python production solver.
6. Paper comparison requires flat Hopkins profile, pressure1e-6, documented softened gravity, effective256² (320² over5²), common density colormap, matching orbit counts. Current baseline is effective170.7² in the paper's4-wide reference region. Paper Fig7 t600 is a milestone; text reports~1000 orbits. Do not silently reduce final goal to t10 or t600. Boundary conditions also require audit: current KH module uses periodic padding, which can create a rotating-background seam in a small matched box; must avoid repeating old boundary-driven failures.

## Questions (one consolidated response requested)

A. Which statements are proven, which are hypotheses, and which proposed remedies would be physically misleading or violate the user's constraints?
B. What is the highest-value next batch (maximum3–4 controlled experiments) using the native CUDA+MPI code and existing checkpoints? Give falsifiable observables, distinguish bug fix from model change, and minimize costly long repeats.
C. Can spin genuinely cure this mechanism at Re∞, or is a coupled deformation/moment or GCL-consistent volume/mass contract required? Connect to the already implemented full-moment prototype only where valid.
D. What may we honestly claim now about long-time stability, and what exactly should be matched to GIZMO before a long run?

Please inspect the small audit JSONs/diffs here as needed and relevant narrow source sections, not the whole project. Return a critical prioritized recommendation in Korean. No code edits or job submissions.
