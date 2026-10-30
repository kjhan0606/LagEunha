# Kepler vacuum free-surface candidate

Date: 2026-10-30
Branch: `codex/kepler-vacuum-free-surface-20261030`
Base: `745e7391cbf8388ec6f905ffb9ba1bbe2b18a38c`

## Scope

This isolated candidate records the sharp annular benchmark contract and
develops an exact gas-to-vacuum Riemann helper. It does not alter the bounded
W=0 reference branch or job `442401`. The helper is not connected to the
tessellator, force kernel, or active Kepler initializer.

## Source findings

- `Exam/mpirks.exam2d.c::paddingTreeVorork4Particles` adds only MPI halo copies.
  It does not construct Kepler inner- or outer-radius boundary sites.
- `Exam/exam_gpu.cu::voronoi_tessellate_kernel` clips real cells by the optional
  fixed polygon. Such support edges are tagged separately from generator faces.
- `getAccVoro2DBlend_kernel` handles ordinary `MAX_INDEX` mirror ghosts before
  the `av_mode == 5` HLLC branch. Its current ghost branch uses the wall-style
  centered pressure/Monaghan path. Under `GFS_VARPRESS_W0`, a ghost-pair face
  aborts as unsupported.
- The HLLC path floors both reconstructed pressure and density at `1e-10`, so a
  literal zero-pressure/zero-density state cannot be passed through unchanged.
- `ParticleSoA` currently carries particle IDs but no boundary-kind array.
  `is_ghost` is derived from the single `MAX_INDEX` sentinel; a vacuum mirror
  therefore needs explicit boundary-kind metadata rather than overloading the
  reflecting-wall behavior.
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
trace; cell-center values are only first-order. This is a local W=0 closure
only. Curvature, multiple faces sharing one ghost site, and nonzero Laguerre
weight rates need additional treatment.

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

## Verification performed

`Exam/tests/test_gfs_vacuum_riemann.c` checks the fan head, the analytic
`gamma=5/3` midpoint `(c,rho,p)=(0.5,0.125,0.01875)`, entropy preservation,
the vacuum edge and zero ALE flux, a shifted `gamma=1.4` state, and invalid
input rejection, including NaN and infinity rejection. It now also checks the
Simpson-integrated fan mass and the nonzero mass, momentum, and energy ALE flux
inside the fan, and the outer- and inner-normal mirror kinematics. Strict C99
compilation and execution pass. A standalone CUDA kernel including the header
kernel including the earlier header revision was reported to compile for
`sm_90`; this was compile-only, not a GPU runtime result. The current login
environment has no `nvcc` or CUDA module, so the new mirror-kinematics helper
has only strict C99 CPU verification. Its test checks an oblique,
unnormalized normal. An independent Python reference checks the same fan,
Galilean shifting, zero front flux, and mass integral; its recorded result is
in `../../../vacuum_face_reference_20261030/reference_result.json`. These
verify local analytic fluxes only; no CUDA runtime integration, mesh motion,
conservation in the finite-volume solver, or Kepler stability is established.
No sharp-annulus initializer has passed an MPI or runtime test.

## Integration gates

1. Add an explicit boundary-kind field to the host-to-device particle/face
   path. Keep MPI halo copies and reflecting walls on their existing behavior.
2. Generate distinct inner and outer mirror sites from the current real
   boundary faces at each RK stage. The planar W=0 helper supplies the target
   normal edge speed and mirror velocity; ensure each ghost site serves only a
   compatible face. Verify the resulting face speed and geometric
   conservation law instead of assuming the reflection gives it.
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
