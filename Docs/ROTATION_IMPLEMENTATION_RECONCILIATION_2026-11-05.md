# Rotation implementation history reconciled with material-cell request

The user correctly recalled a previous internal-rotation implementation. The current baseline's lack of an active spin state is not evidence that no such implementation ever existed. No old code is removed or copied into production by this audit.

## Three distinct mechanisms

1. Geometric face rotation: existing CPU and CUDA Voronoi centroid-offset face-speed term; Laguerre extension via GFS_LAGUERRE_ROTATION. This is not an evolved internal spin state. The previous scalar-weight material-face incompatibility analysis already includes this geometric rotation term.
2. Evolved affine internal motion: `/gpfs/kjhan/LagForce/Kepler_AffineSpin_20261028` has `KeplerAffineStateLagrangian2D`, centroid offsets, three covariance entries, four position-velocity moment entries, total energy, and four RK4 stages. Spin is Bxy-Byx, derived from the evolved velocity moments; affine kinetic energy includes both rotation and strain. CUDA reconstructs face velocities from this state and accumulates traction moments and pressure work. Host RK4 stages and GPU storage/transfer are implemented. This is more than a planned model or postprocessing diagnostic.
3. The selected current production baseline and current experimental acoustic FV candidate do not enable/use that affine state. They must not be described as spin-enabled trajectories.

## What was tested

The retained component audit at `/gpfs/kjhan/LagForce/Codex_moment_spin_backends_20261018/analysis.json` records 1590 valid cells, CPU1/OpenMP4 bitwise equality and CUDA/CPU maximum difference 7.98e-13. Its scope is spin-ledger components, not full trajectory validation.

The isolated Q1 384² H200 startup under `Kepler_AffineSpin_20261028/startup_q1_affine_t001_20261028/log` actually executes the affine path. It terminates with an invalid thermal residual at RK4 stage2. The gravity-quadrupole variant's saved log also terminates at stage2. The target0.001 is not an achieved endpoint. These failures establish that this integrated branch was tried, not that spin physics itself is invalid. The earlier 420² debug runs had wrong-IC issues and are not used here as evidence.

Separate later material-moment work also reports resolved rotation and very short accepted full steps, but no t6/8/10 qualification. The Oct22 handover reports step4 at8.2355880701376864e-5 and442.17s for the next whole step; this is historical handover evidence, not independently revalidated performance in the present audit.

## Consequence for current design

Reuse and inspect the existing affine pressure/torque/work algebra before proposing another spin model. Assess its initialization, bulk/affine/thermal split, generator versus mass-centroid motion, and compatibility between the moment shape and actual shared geometric face. Internal rotation in a velocity reconstruction is not automatically an independent change of geometric face normal. Do not claim the scalar-weight obstruction disappears or that the old spin model solves material boundary compatibility without checking that coupling.

The current user requirement is fixed mass and zero mass transfer through each face. The physical12 net-only fixed-mass branch is comparison evidence, not compliance with that clarified condition. Existing successes remain preserved.

Hashes of inspected sources and logs are stored in the material_lagrangian_design audit directory, `rotation_history_sources.json`.


## Follow-up: native material-work control and archived stage attribution

The unchanged old `kepler_affine2d.h` was called in a new small C audit, compiled with warnings treated as errors. This is not a new simulation or a replacement solver. For a square material cell with uniform pressure, translation, linear shear, and an initially rigid rotating velocity field all give zero thermal rate. Uniform dilation gives Udot=-P V div(u), with absolute identity error2.12e-22. Spin-rate is zero for all uniform-pressure controls. The initial rigid rotation is an instantaneous free-flow control, not a claim of sustained rigid rotation without centripetal forces.

A separate exact uniform-acceleration control isolates a time-integration issue. For unchanged mass and thermal energy, an Euler substage evolves velocity as v+h*a and total energy as E+h*m*v.a. Its recovered thermal energy is U-0.5*m*h²*|a|². RK4's first midpoint predictor has precisely this Euler form. The control with U=1.5e-6, a=1000 and h=5e-6 returns U_stage=-1.1e-5, while the exact compatible energy increment preserves U to roundoff. No floor, reset or new physical variable is needed for that analytic exact kick. This is a mechanism demonstrating that stage negativity alone does not disprove the continuous rotation model; it does NOT yet attribute the archived Kepler failure to this mechanism.

The archived failure diagnostic itself mixes stages: its E0bulk and E0th are computed from the OLD affine state but CURRENT, already-updated bulk velocities. Consequently the printed negative E0th is not evidence that the initial thermal state was negative. The sole saved initial snapshot has zero initialized affine records (written before affine initialization), so it cannot independently reconstruct the full failed trial.

Using the initialization identity and unchanged logged ie, the first archived failed trial separates as:

- starting thermal3.73424874e-10;
- total-energy change-3.35080628e-10;
- bulk kinetic change+8.07952046e-10;
- affine kinetic change+2.94593235e-10;
- resulting thermal-1.06420103e-9, matching the logged failed value.

The gravity-quadrupole variant also fails; it is not a successful repair. Exact pressure-versus-gravity and spatial-versus-stage attribution requires a correctly staged state/rate capture. Preserve these outcomes rather than rerunning the full Kepler campaign blindly.

The user's no-face-mass-transfer direction remains authoritative. The net-only physical12 branch is not selected for continued development. Its complete-mass-flux EC/acoustic energy expression is not automatically the pure material pressure-work expression simply because net mass flux is constrained. Reuse the old affine material pressure/traction/work helper as the local algebraic baseline while auditing actual geometric boundary motion.

Relevant external context: Maire et al. (2007), https://doi.org/10.1137/050633019, couples nodal mesh velocities and pressure fluxes for conservative Lagrangian hydrodynamics. This supports testing force/work/geometry jointly; it neither proves compatibility with a scalar Laguerre weight representation nor justifies adopting its complete method without cost/shear checks.
