# Kepler vacuum free-surface candidate

Date: 2026-10-30
Branch: `codex/kepler-vacuum-free-surface-20261030`
Base: `745e7391cbf8388ec6f905ffb9ba1bbe2b18a38c`

## Scope

This isolated candidate records the sharp annular benchmark contract and
develops an exact gas-to-vacuum Riemann helper. It does not alter the bounded
W=0 reference branch or job `442401`. The working tree now also has
host/device boundary-kind plumbing, but it still has no vacuum-generator
construction or free-surface flux evolution.

## Source findings

- `Exam/mpirks.exam2d.c::paddingTreeVorork4Particles` adds only MPI halo copies.
  It does not construct Kepler inner- or outer-radius boundary sites.
- `Exam/exam.c::mkLinkedList2D_oExam` bins local particles and then MPI
  padding; it does not add physical-boundary mirrors. Vacuum sites therefore
  need a Kepler-specific stage hook, separate from generic MPI padding.
- `Exam/exam_gpu.cu::voronoi_tessellate_kernel` clips real cells by the optional
  fixed polygon. Such support edges are tagged separately from generator faces.
- The current working-tree patch adds `ParticleSoA.kind` and
  `FaceCSR.boundary_kind`. Both fused and ordinary CPU face extraction classify
  a ghost carrying `Wallflag` as a vacuum generator; GPU tessellation maps that
  particle kind to a face kind and handles fixed-polygon sentinels separately.
  The direct CPU force path, GPU force kernel, and CPU reference fallback now
  fail closed on a vacuum face. Host-side integration is syntax-checked, but
  the CUDA kernel has not been compiled.
- Fixed polygon supports arrive at `getAccVoro2DBlend_kernel` as
  `neighbor_idx <= -1000`. They use zero laboratory-frame face speed and a
  wall-work correction, so they cannot represent an expanding vacuum edge.
- `getAccVoro2DBlend_kernel` handles ordinary `MAX_INDEX` mirror ghosts before
  the `av_mode == 5` HLLC branch. Its current ghost branch uses the wall-style
  centered pressure/Monaghan path. Under `GFS_VARPRESS_W0`, a ghost-pair face
  aborts as unsupported.
- The HLLC path floors both reconstructed pressure and density at `1e-10`, so a
  literal zero-pressure/zero-density state cannot be passed through unchanged.
- The active `EUNHA_IC=kepler` path is `Exam/KH/util.c::kh_mkinitial`, not
  `Exam/Kepler/util.c::kp_mkinitial`. It assigns the positive Hopkins floor
  and tapers, then constructs mass and internal energy from that profile. A
  first annulus prototype was mistakenly attached to the legacy Kepler driver
  and has been removed; no sharp-annulus mode is currently implemented.

## Exact local gas-to-vacuum state

For a gas state `(rho0, p0, u0)` on the left and vacuum on the right, let
`c0 = sqrt(gamma p0/rho0)`. The right-going rarefaction head and vacuum edge
move at

```text
xi_head = u0 - c0
xi_vac  = u0 + 2 c0 / (gamma - 1)
```

Inside the fan, at `xi = x/t`,

```text
c   = (gamma - 1)/(gamma + 1) * (xi_vac - xi)
u   = xi + c
rho = rho0 * (c/c0)^(2/(gamma - 1))
p   = p0 * (c/c0)^(2 gamma/(gamma - 1))
```

The new header `Exam/gfs_vacuum_riemann.h` evaluates this state and its ALE
Euler flux. In the face-normal coordinate, use the same equations at the
normal face speed. At the vacuum edge, the mass, momentum, and energy fluxes all
vanish. For an inner annulus boundary, orient the normal from gas into the hole
before applying the same right-going solution. For an equal-weight Voronoi
mirror site across a locally planar vacuum edge, the face speed is
`w=(u_n+u_ghost,n)/2`. Therefore use `w_edge=u_n+2c/(gamma-1)` and
`u_ghost,n=2w_edge-u_n`; the helper returns this edge speed and ghost velocity,
copying the tangential component. Its primitive input must be a gas-side face
trace; cell-center values are only first-order. A second helper solves two
independent face-normal constraints for one shared ghost site, allowing each
adjacent fluid cell to supply a different state. Singular parallel constraints
fail closed. Three or more incident faces, curvature, and nonzero Laguerre
weight rates still need additional treatment.

## GIZMO particle-method initial-condition contract

Reference: Hopkins (2015), §4.2.4 and footnote 24,
<https://arxiv.org/abs/1409.7395>.

The matched particle-method benchmark contract is the cell-centered `256 x
256` grid over a 4-by-4 square centered on the simulation box, retaining the
sites with `0.5 <= r <= 2`. The reference realization has 48,240 sites. Each
site starts with area `1/4096`, surface density `1`, mass `1/4096`, pressure
`1e-6`, and circular velocity in the configured softened Kepler potential.
This is a target specification only; the active initializer still uses the
bounded tapered profile. Applying this sharp cut to the current bounded
tessellator would assign fluid volumes to the vacuum region, so it must not be
enabled before the free-surface geometry and flux closure pass their gates.

### Lattice-completed boundary geometry audit

An independent combinatorial audit of the target lattice gives a low-overhead
starting representation for the vacuum geometry. Keep the 48,240 retained
sites as fluid sites and add the unique cardinal-neighbor lattice sites that
were removed by the annulus cut as *vacuum geometry generators*, not fluid
cells. This yields 908 such generators and 1,280 fluid-vacuum faces (1.88%
generator overhead). Of the 908 generators, 536 touch one fluid site and 372
touch two; in the latter case the two initial face normals are orthogonal.
Thus the regular-lattice Voronoi cells remain squares of area `1/4096` at
initialization, and a shared generator can satisfy both initial normal-face
speed constraints in 2D. This is a geometry count only: it does not validate
the moving-boundary treatment or later topology changes.

The gas lattice extends to within half a cell of `x=±2` and `y=±2`. Eighty-
eight of the missing cardinal neighbors lie outside the nominal `[-2,2]^2`
box, so the computational/tessellation domain needs at least a one-cell guard
margin while retaining the paper's original lattice coordinates. The sum of
the retained square-cell areas is `11.77734375`; the continuum annulus area is
`11.780972450962`, a `-3.08e-4` fractional difference from center-based
particle pruning. The square-cell domain is therefore a resolution-scale
staircase realization, not an exact circular cut-cell geometry.

For the regular target, one guard cell on each side can preserve the physical
spacing and mass exactly: use `NX=NY=258`, a computational box width
`4.03125`, and unchanged 256² target coordinates with spacing `4/256=0.015625`.
Then the initializer's `meanvol=Lx*Ly/(NX*NY)` remains `1/4096`, matching each
fluid cell's square area. Set the Kepler center explicitly to the center of
the original 4-by-4 target lattice. The active initializer's default
`EUNHA_KEPLER_CX=0.5*Lx` would shift the center by half the guard width when
`Lx=4.03125`; that default must not be used for this guarded setup. This is
only a startup-geometry minimum, not a long-run box choice. The independent
planar vacuum-fan reference places the outer front near `r=4.31` at `t=600`.
That no-gravity planar estimate omits Kepler acceleration, rotation, and radial
curvature, so it cannot set the two-dimensional disk's required domain size; it
motivates a domain-sensitivity study. A width-10 box with `NX=NY=640` would
preserve `dx=1/64` and mean area `1/4096` while providing an `r=5` margin. This
is a candidate for that study, not a qualified domain or physical boundary
condition. Disable the active fixed polygon and check periodic-image isolation
and domain convergence before any long run.

The current change is metadata plumbing and a fail-closed guard, not a free
surface. The sparse glass loader still treats each loaded coordinate as a fluid
site; no Kepler stage hook constructs zero-mass vacuum generators. The exact
Riemann helper is not called by the force kernel, and no per-stage generator
motion, gas-side face reconstruction, ALE flux, or cell-volume geometric
conservation check is implemented. A shared generator may only be used while
its incident-face kinematics are mutually satisfiable; the tested pair helper
covers two independent planar constraints, not later topology changes,
curvature, or nonzero Laguerre weight rates. The current GPU path emits NaN on a
vacuum face and the CPU path aborts, so it cannot silently treat vacuum as a
reflecting wall. No MPI/GPU run has validated the boundary representation.

### Verified RK-stage lifecycle

`kh_mkinitial` calls `kh_findVol` to measure the initial volumes. That is only
the initial-condition path. During evolution, `exam2d_vph_rk4_int_blend`
recomputes density/pressure and then force at each RK stage. The update path
calls `mkLinkedList2D` before its cell-volume pass; the force path builds the
linked list again before extracting faces. Both calls use the same current
fluid stage state, but they currently have separate geometry builds.

Consequences for the next implementation:

- Do not create vacuum sites only in `kh_mkinitial` or `kh_findVol`; RK
  predictor stages would then use stale geometry.
- Construct virtual sites from each current RK stage after particle migration
  and MPI padding, and use identical coordinates, IDs, and weights in the
  volume and force builds for that stage. Otherwise the volume used in the
  state update can disagree with the faces used in the flux, invalidating the
  geometric conservation check.
- Keep virtual sites outside `VORO_NP`: that count drives fluid updates,
  conserved output, and migration. A stage-local buffer must distinguish
  vacuum generators from ordinary MPI halo copies and be released before the
  next stage rebuild.
- The linked-list callback is a convenient insertion point for already-known
  virtual sites, but it cannot infer a deformed annular boundary from a static
  lattice mask. Boundary incidence and site motion still need a stage-aware
  construction rule; the 908-site stencil audit proves only the initial grid
  topology.

The reproducible count and guard-width audit is
`Exam/tests/audit_kepler_annulus_vacuum_stencil.py`.

## Verification performed

`Exam/tests/test_gfs_boundary_kind.c` passes strict C99 compilation and
execution. `mpiicx -fsyntax-only` passes for `exam.c` and
`exam_gpu_extract.c` with `USE_CUDA` enabled, covering the host-side paths; it
reports two unrelated existing warnings in `exam.c`.

`Exam/tests/test_gfs_vacuum_riemann.c` checks the fan head, the analytic
`gamma=5/3` midpoint `(c,rho,p)=(0.5,0.125,0.01875)`, entropy preservation,
the vacuum edge and zero ALE flux, a shifted `gamma=1.4` state, and invalid
input rejection, including NaN and infinity rejection. It now also checks the
Simpson-integrated fan mass and the nonzero mass, momentum, and energy ALE flux
inside the fan, plus outer, inner, and oblique-normal mirror kinematics. It
also solves a shared orthogonal generator from two different gas states and
confirms both edge speeds give zero ALE flux; parallel constraints are
rejected. The stencil audit checks the initial lattice counts above. Strict C99
compilation and execution pass; the derived outer and inner edge speeds also
yield zero ALE mass, momentum, and energy flux. A standalone CUDA kernel
including an earlier header revision was reported to compile for `sm_90`; this
was compile-only, not a GPU runtime result. The current login environment has
no `nvcc` or CUDA module, and its Slurm node inventory exposes no GPU GRES, so
the current CUDA source changes have not been compiled or run. The new
mirror-kinematics helper has only strict C99 CPU verification. Its test checks an oblique,
unnormalized normal. An independent Python reference checks the same fan,
Galilean shifting, zero front flux, and mass integral; its recorded result is
in `../../../vacuum_face_reference_20261030/reference_result.json`. These
verify local analytic fluxes only; no CUDA runtime integration, mesh motion,
conservation in the finite-volume solver, or Kepler stability is established.
No sharp-annulus initializer has passed an MPI or runtime test.

## Integration gates

1. Add a Kepler-specific RK-stage hook that generates inner and outer
   zero-mass vacuum sites from current real boundary faces. Keep MPI halo
   copies and reflecting walls on their existing behavior.
2. Generate inner and outer mirror sites from current real boundary faces at
   each RK stage. The initial Cartesian mask needs one site per missing lattice
   coordinate: shared sites each have two orthogonal face constraints, which
   the pair helper can satisfy from separate gas traces. Rebuild or reclassify
   this incidence after the mesh distorts; fail visibly for incompatible or
   over-subscribed sites. Verify actual face speeds and the geometric
   conservation law instead of assuming the mirror construction gives them.
3. Add a vacuum-face flux branch that bypasses positive HLLC floors and applies
   the exact gas-to-vacuum state. Keep real-cell mass fixed and account for
   moving-face pressure work consistently.
4. Qualify a planar expansion into vacuum first: front position, rarefaction
   profile, total energy, and face-volume/GCL residual. Then use a coarse
   annular Kepler test with both boundaries before attempting 256^2 or a long
   continuation.

The smallest missing solver contract is boundary-face geometry and flux
consistency together. A vacuum state alone or a pressure floor alone does not
create a free boundary.
