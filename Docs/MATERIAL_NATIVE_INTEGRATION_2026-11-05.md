# Native 2D material-work integration: first implementation stage

The existing `Exam/exam.c`, `Exam/exam_gpu_extract.c` and `Exam/exam_gpu.cu` now call `gfs_face_work2d` from a shared C/CUDA header for the ordinary total-energy traction work. Existing velocity selections and pressure models are unchanged. The variational Kepler pressure branch is unchanged. This is an implementation foundation, not a new material-boundary evolution scheme or a physics fix.

The helper exposes oriented geometric sweep and work-velocity versus mesh-velocity power difference. The difference is diagnostic only; it is never added as an energy correction. `Exam/Makefile` includes the new header dependencies.

Validation performed on the Ada host:
- Actual native MPI/OpenMP C source, extraction C source and CUDA source compile independently using the repository Makefile flags. All three exit 0; existing warnings remain. Objects are isolated in `/gpfs/kjhan/LagForce/material_native_integration_20261105_v1`, preserving existing binaries.
- Native C convex-polygon test compares shoelace area differentiation against integrated edge motion: discrepancy 3.44e-12; uniform-pressure work discrepancy 7.91e-12. Shared-face traction work cancels pairwise. Unequal pressures expose a nonzero work defect despite cancelling sweep defects.
- Actual Ada kernel evaluates 4096 faces; maximum normalized host/device work and defect discrepancy 5.3291e-15 with the production fast-math CUDA flag.

No new full 2D hydrodynamic time evolution has been launched. MPI distributed evolution, mesh-motion compatibility, positivity/convexity through evolution and KH/RT/shock accuracy remain unverified for the proposed material update. Existing `GFS_GEOM_FACE_VEL` already selects geometric work velocity; merely activating that old option would not implement the successful 1D contact-speed material boundary in 2D. Next integration must connect shared mesh motion and work, rather than just choose a different work velocity or introduce mass exchange.

## Shared-vertex kinematics control

Added `Exam/gfs_material_mesh2d.h` with host/device area and exact straight-line edge sweep, indexed by shared vertex IDs. The strict-convexity check examines the quadratic edge/vertex orientation throughout each proposed linear-motion step, including interior minima; it rejects degeneracy and inversion rather than clipping or steering. The caller must reject the complete mesh step transactionally.

The native C geometry control performs 1000 prescribed-shear steps to t=1 on two adjacent convex cells. Shared-edge sweeps cancel, fixed-mass densities remain constant, and area changes agree with edge sweeps to 1.249e-16 including an additional compressing/rotating affine step. Deliberate inversion is rejected. This is NOT a hydrodynamic time evolution: prescribed vertex velocities are inputs, pressure does not determine those velocities, and existing FaceCSR ownership/MPI integration remains pending. No KH/RT run is claimed.

The NativeU64 production-source counterparts have also been connected to the common traction-work helper and compiled successfully as C/CUDA objects with their retained feature flags. Build commands/results are in `native_build_results.json` in the isolated integration directory.

## Pressure-driven shared-node operator and coupled evolution controls

Implemented the inverse of the retained exact ideal-gas wave curve in `gfs_exact_riemann_ale.h` as a host/device piston-pressure function. `gfs_material_node2d.h` solves the two-component shared-vertex pressure-force balance over incident half edges. Every cell contribution uses that same solved vertex velocity for traction work. No fitted viscosity coefficient, mass flux, thermal reset, or independent velocity smoothing is introduced. A safeguarded local Newton solve reports failure rather than silently falling back.

Four planar reductions (Sod, strong pressure jump, double rarefaction, colliding shocks) reproduce the original exact Riemann star pressure/velocity within 9e-16 normalized differences. An asymmetric four-cell vertex produces a nonzero pressure-driven velocity, force residual 7.02e-16 and power residual 4.44e-16. Irreversible heat relative to cell-mean pressure work is nonnegative in these checks. A Galilean shift test passes. An actual Ada kernel solves 256 boosted nodes with maximum velocity discrepancy 5.55e-16.

`Exam/HydroCore/material_node2d_evolution.c` is an integration-test fixture for these operators, NOT a replacement production simulator. It uses shared quadrilateral connectivity, immutable cell masses, cell-mean states (FIRST ORDER IN SPACE), reflecting walls and SSPRK2. Pressure determines vertex velocities; positions, momenta and energy actually evolve together. Every stage checks positive area/pressure and strict convexity. There are no geometric remaps. The new production CSR/MPI ownership integration and second-order reconstruction are still absent.

Actual coupled time evolution:
- 8x8 static contact and smooth pressure pulse both reach t=0.1 in 59 steps, no rejected steps. Contact remains stationary with uniform pressure. Pulse energy error 3.51e-16.
- 16x16 counterparts reach t=0.1 in 118 steps, no rejected steps. Contact energy/area/velocity errors are zero in the measured output. Pulse energy error -7.02e-16, area error 4.44e-16, minimum pressure 1.00011869. Maximum node balance residual 1.99e-13. Recorded native CPU runtimes approximately .05/.10 seconds for these tiny controls only; no production performance inference.
- Evidence: `/gpfs/kjhan/LagForce/material_native_integration_20261105_v1/node_evolution_controls_v1/results.json` and per-case stdout/stderr; full native-source builds and prior kernel tests remain in the parent directory.

These checks establish a pressure-driven coupled 2D control, not KH/RT/Sedov validation, not a second-order method, not a production MPI/GPU trajectory, and not equivalence of material polygons to instantaneous power cells. Smooth shear dissipation, strong-shock positivity, topology evolution and distributed ownership remain unqualified. No existing simulation was restarted or replaced.

## Actual retained Voronoi/Laguerre geometry: import and failed evolution gate

Added a bounded, single-rank FaceCSR/ParticleSoA bridge in the existing ordinary and NativeU64 extraction source. Roundoff-scale spatial hashing identifies shared vertices; ambiguous matches, unmatched interior edges, inconsistent loops and non-convex polygons fail explicitly. MPI/periodic images are explicitly unsupported pending ownership work. The import is not enabled as a production evolution mode. No new tessellator was written.

The test calls the retained `gfs_fv_geometry.c` / `Voro/voro.c`, imports their real faces, and checks native cell areas. At 64 cells, regular Voronoi has 81 shared nodes; slightly displaced Voronoi and weighted Laguerre have 130. All fill unit area; maximum native/imported cell area discrepancy 1.39e-17; largest merged coordinate distance 1.57e-16.

Pressure-driven evolution on these actual imported polygons reveals an unresolved convexity failure that the rectangular controls missed:
- Voronoi, constant reconstruction: stop t=0.03182853865; dt/2 stops t=0.03183039072.
- Weighted Laguerre initial mesh, constant reconstruction: stop t=0.01039405298.
- Limited linear reconstruction: Voronoi stop t=0.04514461721; weighted Laguerre stop t=0.04634152287.
- All tested stationary contacts reach t=0.1, with uniform pressure and negligible velocity.

At the first Voronoi stop, total area is 0.9999999999999999, minimum pressure is 1.00000478 and the problematic cell has positive area 0.01484296. One orientation determinant is 1.30e-17 with derivative -4.16e-4: the cell is about to become concave. The nearby short edge still has finite length 2.17e-5; a zero-length-edge reconnection is not already available. Step halving approaches the convexity boundary rather than removing it. Stored states and `native_pulse_failure_analysis.json` provide the evidence. The last accepted state is preserved; no clipping, steering, remap, mass transfer or convexity relaxation was used.

The linear reconstruction uses centroid least squares, scalar vertex extrema bounds and the retained `gfs_velocity_midpoint_alpha` limiter. Observed second-order convergence has NOT been established. Its positive trace-level acoustic dissipation does not prove cell-mean entropy stability. Therefore neither reconstruction order qualifies for KH/RT yet.

The user explicitly authorized completion through BOTH KH and RT without individual approval. Those tests remain pending the convexity/representation gate, not pending user permission.

## Why exact material boundaries do not guarantee convexity

An elementary inviscid Euler control clarifies the distinction. On an unbounded fluid, take constant density/pressure and stationary velocity u=(0,x(1-x)). Divergence and material acceleration both vanish. The exact material map is x=X, y=Y+t X(1-X), with Jacobian one. For the initial unit square, the two lower corners remain (0,0) and (1,0), but their midpoint (.5,0) lies below the transported lower boundary y=t/4. Thus this material set is non-convex for any positive time despite exact mass conservation and zero heating. This is an analytic statement about exact material sets, not a proof that approximate fixed-mass Laguerre hydrodynamics is impossible.

A globally convex numerical tessellation is consequently an approximation to material geometry, not an automatic consequence of zero mass exchange. The current freely advected vertex prototype cannot be adopted under the user's mandatory convexity requirement. A return to the retained convex Laguerre representation must be distinguished from exact material-boundary tracking and must retain consistent discrete volume/pressure work; merely retessellating these saved states and silently reassigning their contents is not a solution.
