# Fixed-mass shared-face evolution: first implementation

The user requested solving the incompatibility, not stopping at repeated diagnosis. An isolated native C 1D material-cell driver now reuses the existing exact Riemann solver; shared face contact speed updates both boundary positions and pressure work. Momentum uses the same shared pressure. Mass is immutable. No new stress, sensor, floor, or heat repair is introduced. Spatial first-order Godunov plus SSPRK2. This is a comparison implementation, not a replacement for the user's native OpenMP/MPI/CUDA production code.

Eight runs (N200/N400) reach their full registered endpoints: moving contact t=.25; Sod t=.2; double rarefaction t=.15; reflecting Blast t=.038. Every cell mass stays bitwise identical and neighbor faces remain exactly shared. Energy balance including boundary work is within 1.4e-14 relative to max(1, initial energy); positive cell widths/thermal energy retained.

Moving-cell-average density L1 errors (N200 -> N400): Sod 2.345% -> 1.472%; double rarefaction 1.839% -> 1.275%; Blast 16.178% -> 10.510%. Moving contact remains at roundoff. Sod/doublefan reference is the analytic Riemann fan sampled on 20000 bins; Blast uses existing independent Athena N6400 finite-resolution reference. Blast density spikes and diffusion remain; completion and conservation do not constitute full accuracy success.

This establishes a consistent 1D material-boundary path, not a solution to multidimensional Laguerre geometry. Noh, Shu-Osher800, 2D/3D convex shared material geometry, long Kepler, rotation coupling, and matched GIZMO performance remain open. The exact iterative solver and tiny 1D timing do not qualify production speed. Prior weighted/zero-weight results and all checkpoints are preserved.

Evidence: /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_shared_face_evolution_20261105_v1

All eight simulations completed. No ongoing simulation is claimed.

## Second-order reconstruction comparison

Evidence directory: /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_shared_face_plm_20261105_v1

Parameter-free minmod primitive PLM on nonuniform cells, same shared contact velocity/pressure work and SSPRK2. Eight first-order regressions have byte-identical snapshots. All eight PLM endpoints pass. Moving-contact preservation, per-cell fixed mass, shared face coordinates and energy conservation remain intact (energy error <1.6e-14).

At N200 density average L1: Sod 2.345% -> .689%; doublefan 1.839% -> 1.019%; Blast16.178% ->4.998%. PLM N400 errors .357%, .701%,2.554%. Smooth acoustic N64/128/256 gives observed rates1.876/1.903, compared with first-order .987/.993, against a small-amplitude linear wave reference.

Do not equate improved global L1 with eliminating artifacts: Blast N200 density max increases8.035 ->8.152; PLM N400 max7.250. Sod contact-window max pressure error at N200 increases slightly .01644% ->.01826%, while PLM N400 falls to .00769%. The narrow density spike remains unresolved. No extra damping patch added. Multidimensional shear, convex material geometry and native long Kepler remain unverified.

## Resolution-sensitivity follow-up

User observed nearly overlapping second-order curves. All sixteen added N100/N800 runs (both orders, four cases) completed; existing N200/N400 runs reused unchanged. Independent interval overlap measures field changes without relying on plotting interpolation or a reference solution.

Relative density-field change when doubling200->400: Sod first-order0.9311%, second-order0.4584%; double rarefaction1.5216% vs1.3505%; Blast6.6972% vs3.0759%. Thus second order is substantially less sensitive for Sod/Blast, but the reduction for double rarefaction is modest. No universal resolution-independence claim.

Second-order density mean L1 error N100/200/400/800: Sod1.2856/.6883/.3564/.1833%; doublefan1.7923/1.0186/.7007/.4761%; Blast10.4587/4.9980/2.5543/1.7711%. The error still depends on N and decreases on refinement. Sod and Blast second-order N100 roughly match or improve first-order N400 global density error. This is accuracy equivalence for these fixtures, not a runtime claim.

Important exception: Blast density peak7.4176/8.1516/7.2497/7.8228 is nonmonotonic with N. Narrow peak artifacts remain unresolved despite converging bulk profiles. Contact is preserved at roundoff for both orders and thus does not distinguish their accuracy.

Analytic reference sampled on40000bins (previous20000), integrated on each actual cell; small differences from prior quoted errors are quadrature resolution. Blast reference is Athena N6400, not exact; its previous N3200/6400 density discrepancy is1.1267%, so the1.7711% comparison atN800 should not be interpreted as precision against an exact solution. Adjacent-simulation field differences are independent of this reference uncertainty.

All cells retain original mass bits, shared boundaries and positive widths. This study remains1D and does not establish cold shear, rotating-cell or long-Kepler robustness.

Evidence: /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_shared_face_plm_20261105_v1/resolution_extension
