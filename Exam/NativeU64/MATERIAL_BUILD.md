# Native material build profile

The opt-in `material` profile connects the common periodic, reflected-wall,
and prescribed-exterior geometry in the existing native C/CUDA code.
The default profile remains `legacy`, including its Kepler geometry settings.
This profile currently requires one MPI rank. It is not a completed 3D or
multi-rank hydrodynamic implementation.

To compile the six material modules from the current checkout:

```sh
make -C Exam material-native-objects \
  MATERIAL_BUILD_DIR=/gpfs/your-project/fresh-build-directory CUDA_ARCH=89
```

Use an absolute, nonexistent output directory and the architecture of the
target GPU (89 for the tested Ada, 90 for H200). Intel MPI compilers, CUDA,
and MKL must be available in the environment. The face-response and weight
kernels retain the validated sm_80 target; architecture retuning is separate.
The material weight solve uses the existing cooperative persistent GPU PCG.
It retains FP64 arithmetic, the projected Jacobi preconditioner and the
independently recomputed relative residual threshold of `1e-10`. The runtime
checks cooperative-launch support; unsupported devices fail explicitly.
History recycling is not enabled by the build profile. Current Ada checks
and late KH/RT comparisons are recorded in `native_persistent_weight_v1`
under the material pressure-acceleration research directory. This backend
selection does not qualify untested hardware, MPI layouts or physics.
The target writes its commands, object hashes, and a module archive. It does
not by itself rebuild the rest of the executable's dependency archives.

For a full source build, use a fresh checkout without existing objects or
archives and a Slurm compute allocation, then run:

```sh
NATIVE_PROFILE=material BUILD_JOBS=2 bash Exam/NativeU64/build.sh
```

The default target is H200 (`CUDA_ARCH=90`). For CPU-only compilation on an
authorized development host, explicitly use `BUILD_LOCAL=1`; for subsequent
execution on Ada, also set `CUDA_ARCH=89`. No GPU is used during compilation.
The script keeps the fresh-checkout requirement. A link-only resume must use
the same profile as the compiled objects. Adding this build route does not
qualify its physical results; source-module integration comparisons are recorded in
`Docs/MATERIAL_RT_PORT_AUDIT_2026-11-06.md`.

All runtime physics options and restart provenance must accompany a result.
A successful compile or short restart comparison does not establish spatial
convergence, long-time Kepler accuracy, or the matched-accuracy GIZMO timing
target.

The material profile also dispatches the density/pressure gradient projection
on CUDA. Host and device include the same scalar reconstruction and
1/2/3-dimensional constrained-gradient solver; no copied device formulas or
new physical coefficients are introduced. This kernel retains FP64 with
`--fmad=false`. The existing prescribed-exterior path is unchanged.
`GFS_CUDA_PROJECT_AUDIT=1` enables an expensive CPU shadow comparison for
qualification, not production timing. In `native_common_projection_v1`,
matched KH/RT/Sedov continuations agree to rounding-scale differences; Noh
is an unchanged-path nonregression check. Two late-KH paired timings show
4.5--4.8% wall reduction on the shared Ada host. These checks do not qualify
3D evolution, multi-rank MPI or the full matched-accuracy performance target.

The periodic weight RHS, its geometric transpose, and the prescribed-exterior
face-rate wrapper now use `Exam/gfs_material_power_face_motion.h`. For a
directed face of measure A, site separation d, and owner-relative face
centroid c, the common geometric volume rate is

```
Vdot_face = (A/|d|) [c.vi + (d-c).vj + (Wi_dot-Wj_dot)/2].
```

The transpose uses these same coefficients; it does not introduce a separate
pressure-force approximation. In 1D, A=1; in 2D it is edge length; in 3D it
is face area. A curved image xj=g(x_owner) requires its coordinate derivative
J in the velocity and J-transpose in the force, rather than the fluid
velocity-reflection matrix. These are geometric relations, not a new mass
flux or a claim that a curved reflected-image wall is stationary.
The dimension-generic C/CUDA algebra is checked separately from the actual
2D native continuations in `common_face_motion_v1`. Three-dimensional algebra
checks do not qualify 3D hydrodynamic evolution.

`gfs_material_face_normal_prepare2d_wall_velocity` accepts an optional array
of normal boundary speeds, indexed by trace face. Only one-sided direct-wall
faces read it. Their ordered-jump limiter operates on velocity relative to
the boundary; the fluid velocities themselves are not modified. Callers must
use the same boundary speed for geometry and moving-wall pressure work.
The existing `gfs_material_face_normal_prepare2d` wrapper supplies NULL and
retains stationary-wall behavior. The nonzero-speed route is currently tested
with the existing C trace builder and moving-wall Riemann function on saved
Cylinder geometry; it is not yet a native Cylinder evolution route. Frozen
trace checks and native stationary/periodic regressions are recorded in
`cylinder_trace_kinematics_v1`. A correct limiter does not resolve a mismatch
between cell-mean fluid velocity and the velocity of a mesh generator.


Moving-wall RHS APIs are now available in both `gfs_material_rhs2d.h` and
`gfs_material_cuda_rhs.h`. Their caller supplies the SAME per-face normal
speed used by `gfs_material_face_normal_prepare2d_wall_velocity`; the kernels
use it for the reflected Riemann state, Vdot and total-energy work and return
an external boundary-power ledger. Legacy stationary entry points remain.
`gfs_material_image_weight2d.h` is an optional anchored image-weight adapter;
it eliminates matching boundary volume terms, then recovers image speeds.
It requires resolved periodic aliases and a prescribed exterior anchor for
every non-image graph component. It is not a stationary cylinder solver.
Random C/CUDA parity, native unchanged-path restart regressions, and frozen
saved-cell annulus coupling are in `cylinder_moving_rhs_v1`. The actual native
Cylinder evolution driver has NOT yet been connected or qualified. A moving
image wall's energy ledger must not be described as stationary wall work.


`gfs_material_image_trace2d.h` and `gfs_material_image_rhs2d.h` now couple
one-sided image walls, prescribed exterior samples, the existing reconstruction,
Riemann volume targets, image-weight solve, and the same CPU/CUDA RHS with
recovered boundary speeds. The optional physical-plane CSR entry point keeps
actual wall geometry, while legacy reflected-wall APIs retain bisector checks.
The adapter returns separate wall/exterior energy ledgers and does not advance
state. It currently performs two RHS evaluations and has no performance
qualification. `cylinder_coupled_trace_v1` tests actual saved-cell annuli at
three grids with equal/nonuniform weights, plus unchanged native continuations.
This remains an instantaneous mixed-boundary implementation: Cylinder native
stage-loop mapping/injection/rollback is unported, and a moving radial-image
wall is not a validated stationary obstacle. Do not bypass qualification gates.

`gfs_material_cylinder_boundary2d.h` now maps the original Cylinder padding
metadata to separate fluid/coordinate velocities and wall/weight/periodic
ownership. C finite-difference and CUDA parity checks are in
`cylinder_original_padding_kinematics_v1`. This adapter does not enable native
Cylinder evolution. Original planar open mirrors have stationary geometric
faces and copied owner weights; they cannot simply be assigned the nonzero
normal speed of a material Riemann face. An unbalanced inflow/outflow exposes
a global volume-compatibility residual even when the geometry is exact. See
the 2026-11-08 boundary audit in `Docs/MATERIAL_RT_PORT_AUDIT_2026-11-06.md`.


`gfs_material_cylinder_rhs2d.h` now adapts the existing mixed periodic,
radial-image and prescribed-exterior operator to native `ParticleSoA`.
Its caller supplies certified alias/image ownership and the same coordinate
low parts used in tessellation. Radial mesh velocity is Jv, while reflected
fluid velocity uses the reflection map. Prescribed exterior m,U,V,v,W remain
boundary data. Both native AoS-to-SoA packers now carry the exterior V64 along
with U64. This does not enable the full RunCylinder route: stage-time, padding
ownership/precision, accepted inflow/departure, wall ledgers and checkpoint
state still need their coherent native connection.

`cylinder_native_bridge_v1` compares the extracted helper against the preceding
operator and compares its C/CUDA evaluations. On the same accepted t=1 Cylinder
state, volume/weight derivative differences are below 5e-16. A t=1-to-1.05
continuation including ownership events preserves exact IDs/masses and agrees
in density to relative L1 ~1.4e-15. This is a short consistency check, not the
400x80 morphology or <=2xGIZMO performance qualification. The native material
objects and complete executable were rebuilt with the current GPUContext.
