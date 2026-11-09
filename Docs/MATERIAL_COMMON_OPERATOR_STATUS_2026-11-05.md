# Common material-cell algebra: current verified scope

The user requires fixed cell mass AND zero material transfer through each face. Energy still transfers through pressure work; zero mass transfer must not be described as zero energy exchange.

Existing 2D affine moment source/traction algebra was generalized by dimensional index bounds to1D/2D/3D. The state contains mass, cell momentum/energy, covariance and momentum-position moments; spin is derived from their antisymmetric part, not a separate heat or spin reservoir. No changed native production source or evolved trajectory is claimed.

Ten analytic uniform-pressure controls pass. Translation and initial shear/rotation have zero thermal rate. Dilation gives dU=-P dV to4.49e-19 absolute. The2D result matches the archived helper to4.34e-19. These tests use exact simple material-cell geometry and do not resolve how native Laguerre faces will follow the chosen material motion.

Fifteen independent normal-Riemann controls reuse the existing exact ideal-gas solver, sampling at its contact velocity. In1D/2D/3D (including oblique normals), the flux reduces to mass0, momentum P*n, energy P*u_normal to3.00e-14. Density contacts remain static at equal pressure. Stationary RH shocks use upstream Mach1.05 and3, with downstream0.9531254 and0.4751910. Reflecting-wall impact at Mach0.3 and3 gives Pwall/Pinitial1.50234976 and17.20832932; stationary-wall mass and energy flux are zero to rounding. This is a face-level reference, not a finite-resolution captured shock or an evolved wall collision.

## Unfinished full objective

- Geometry preference is settled: the user requires material-following general convex cells and mandatory convexity. The unresolved part is a discrete algorithm that realizes this together with fixed mass, conservation and acceptable cost; do not treat the geometric question itself as unanswered.
- Couple pressure/traction, internal moments, geometric volume and a positive compatible time update without adding ad hoc physical damping.
- Re-run original1D endpoint/accuracy tests, including contact-density peaks, with that selected material model. Prior physical11 successes are variable-mass comparisons; physical12 constrains only net mass.
- Demonstrate finite-resolution weak/strong shocks and hard walls in1D/2D/3D; current algebra checks do not establish this.
- Resolve the separately captured3D face-extraction crash; its actual backtrace is in acoustic_residual_energy_dimensional_v1/wall_crash_diagnostic_v1/gdb.log.
- Demonstrate smooth cold shear and rotating flows, resolution/time convergence and relevant longer trajectories.
- Benchmark matched end-to-end runtime against GIZMO and satisfy<=2x. No interface timing or historical component parity proves this requirement.

No additional simulation or agent is running from this bundle. All local native audits terminated. Previous fixed-mass and variable-mass outputs are preserved.

Evidence root: /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_lagrangian_design_20261105_v1

## Native 1D fixed-mass regression added 2026-11-05

Built the current Hydro1DExam/laguerre_sod.c against the common gfs_variational_core.h and ran its nine-problem unified suite at nominal N=200 in an isolated /gpfs directory:

/gpfs/kjhan/LagForce/KeplerVariationalCoreUnify_20261105_v1/oned_regression_20261105_v1/run.log

All nine reached their configured endpoint with positive cell volume and internal energy; the process exited 0. Fixed-mass W=0, global KDK, and the reconstructed-strain alpha=1/beta=2 pair term were used. This is endpoint and positivity evidence only, not an accuracy pass: pressure L1 is 3.56 for Blast and 9.52 for colliding shocks, while Shu-Osher density L1 is 0.866. Stationary and moving contacts remain near machine precision. The test did not measure runtime against GIZMO and does not qualify a multidimensional solver.

H200 CUDA compile/parity job 451241 is still pending with Reason=QOSGrpGRES; Slurm reports a scheduled start of 2026-11-06 22:49 KST.

## 3D face-cycle replay follow-up

The captured 3D Sod wall-failure fixture was rebuilt with bounded face extraction and replayed at N=8^3. The previous memory-corruption crash is now rejected as a geometry error, but the evolution still fails: it reaches t=0.27801238 of 0.32 after 56 accepted steps and 75 retries, then exhausts retries in final geometry. The first logged failure is cell 37, face 8, 21 vertices; the bounded extractor returns an invalid-topology status, not an output-capacity status. Energy-ledger diagnostics remain small, but this is not a 3D pass.

A GDB replay identifies the broken incidence: the face-8 walk reaches a vertex with related planes [20,9,3] from a vertex [9,8,2]. The target does not contain face 8, so the face cycle cannot close. This localizes the defect to Voro3D vertex connectivity, before surface quadrature. A strict audit at trim entry then found an active vertex already linked to an inactive target (mp=22, vertex=5, slot=2); no active replacement shares edge planes (9,8). The dangling link therefore exists before active-vertex compaction, and the old trim routine would index its active-only Old2New map with that inactive target. The caller frame identifies candidate neighbor plane 23 as the cut preceding the failure. A direct reciprocal-link check in new-vertex creation did not trigger, so the exact defect is within the cut/inactivation boundary update. The new guard rejects the geometry before the undefined mapping. No active endpoint exists to relink to, so a heuristic repair is not justified; the clipping/inactivation topology update remains unresolved. The exact-macro replay and audit build are under /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/acoustic_residual_energy_dimensional_v1/bounded_face_linkaudit_residual_20261105_v3; the caller trace is under /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/acoustic_residual_energy_dimensional_v1/bounded_face_cliptrace_20261105_v2. Replay artifacts are under `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/acoustic_residual_energy_dimensional_v1/bounded_face_replay_20261105_v2/replay_sod3d_n8` and the GDB trace used the same input in a fresh GPFS directory.
