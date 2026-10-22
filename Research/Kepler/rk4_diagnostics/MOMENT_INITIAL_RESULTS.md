# Native moment initialization and the next Euler spatial representation

2026-10-20. This bundle is progress toward the active Kepler goal, not a
long-time stability result. All evolution remains native C/CUDA/MPI. The
preceding opinion-only turn did not change implementation or runtime state.

## Implementation completed

Canonical worktree:
`/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`,
origin `git@github.com:kjhan0606/LagEunha.git`, HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`
plus existing and new uncommitted driver edits. The supplied backup cwd is not
this repository. No commit/push or shared-context change was performed.

New code:

- `Exam/gfs_moment.h`: fixed-reference-origin m,q,Q,p,D,E algebra, affine
  kinetic energy, origin/boost transforms, canonical face moment fluxes,
  Euler volume sources and same-point gravity force/torque/work helpers.
- `gfs_moment_density.h`: positive exponential-quadratic density matching
  all m,q,Q moments on the actual polygon. Whitened moments, Newton/Cholesky
  with line search, no moment projection/Hessian regularization/density floor.
- `gfs_moment_initial.h`, `gfs_moment_gpu.inc`, `gfs_moment_native.inc`:
  CUDA integration of the existing continuous IC, reconstruction and
  `GFS_MOMENT_INIT_ONLY=1` gate. It always terminates without evolution.
- The tagged 192-byte payload occupies only the 24-double legacy RK scratch
  prefix, excluding w2backup. Native full-record halo/migration/checkpoints
  keep their 576-byte ABI. Legacy restart/evolution rejects tagged states.
  The gate supports initial W0 geometry only, not a dynamic-weight run.

Nothing in this bundle evolves the proposed moment transport or solves the
weight–mass nonlinear equation. The exact-at-W0 relation is only the initial
identity m=m0; it does not validate a dynamic weight solver. The pressure-only
legacy force path remains separate. A final guard revision adds no default
per-step MPI collective.

## Algebra and actual GPU evidence

`moment_initial_v5_build/contract.json`:

- Six-monomial density realization: pointwise relative error 3.89e-15;
  normalized moment residual 1.33e-15 in four Newton iterations.
- Affine velocity matrix error 1.11e-16; exact volume-source versus direct
  quadrature error 1.39e-16; canonical pair angular residual 2.78e-17.
- Origin/boost invariance, an outside-generator convex cell and an SPD but
  support-unrealizable target are checked. The latter is rejected.
- Actual Euler HLLC massflux versus normal mesh speed is strictly decreasing
  in 1206 samples. Mode bitmask30 covers HLLC/upwind modes1–4, **not HLL
  fallback5**. This is not a Jacobian/SPD proof or a coupled MPI solve.

`moment_initial_v1/analysis.json`: 4096-cell distorted native rigid, shear
and boosted-shear IC integrations all pass. Maximum recovered B error is
2.01e-12 in the boosted case. Thermal subtraction relative error there is
3.17e-8, distinguishable from the actual Kepler representation error below.
Native snapshot payloads match their audits bitwise.

`moment_initial_v5/{manifest,analysis}.json`: all **71080** actual refined
cells pass positive density/mass/S/U and density realization. Maximum
normalized density moment residual is 5.00e-12, at most four Newton steps.
The complete moment payload survives native snapshot writing, reading and
MPI1→MPI4 migration **bitwise**, including the reconstruction cache. Initial
input SHA is unchanged. The unsupported legacy restart is explicitly rejected
with exit99; that is the intended negative check.

Manual user-authorized lageunha/Ada execution: fresh initialization2.756s,
MPI4 restart4.244s, negative legacy guard2.026s. No native calculation/job from
this bundle remains running after completion.

Frozen physics gate binary v5 SHA256:
`61a46b1be350b5e190de2d9635590fe5033461abdee70405fa4ea7ed6277a050`.
`moment_initial_v6_build` contains the subsequent legacy-guard/build-dependency
revision; its initial kernel is unchanged. `moment_legacy_v6/manifest.json`
records the separate unchanged-evolution check.

## Geometry failure, diagnosis and correction

v1–v4 are retained. The actual IC first failed four new positive-triangle
checks, IDs34086–34089. Diagnostic-only reruns were bounded initialization,
not discarded/restarted physical evolution. Four triangles were rejected
because the raw clipping polygon contains a one-coordinate-ulp backwards
micro-edge: x=.046875, dy=-8.67e-19. Existing raw geometry deliberately permits
this inside an explicit arithmetic error bound; an exactly-zero-edge exception
alone did not fix it. No density or thermal failure was concealed.

The new initializer computes vertices from intersections of the same adjacent
W0 power-support planes instead of successive segment interpolation. This
produces65536 exactly inactive edges on the cocircular grid, all of zero
quadrature measure. No edge-length cutoff, weight clipping or physical
parameter was introduced. Every nonzero centroid-fan triangle remains positive.
The initializer's construction is private: a future evolving method must use
the same support vertices for its volume and face quadrature. The old force
path has not silently received different geometry.

## New physical limitation: affine energy is insufficient for this cold IC

All quantities below use the complete successful v5 IC, not the partial v2
sum. Mass75.97420601023128; true integrated thermal sum8.639995880128927e-4.
The same positive quadrature integrates the original curved Kepler velocity,
momentum moments and full total energy. Independently integrating the squared
velocity residual confirms

    E - K_affine = U_physical + K_nonaffine.

The identity error is at most7.90 eps times the kinetic/total subtraction
scale. Thus the extra term is physical velocity representation error, not
roundoff or evidence for adding heat to old depleted checkpoints.

| Initial region | median K_nonaffine/U_physical | maximum |
|---|---:|---:|
| r<eps=.047 | 7.40% | 10.02% |
| .047<r<.1 | .984% | 4.26% |
| .1<r<.2 | .0767% | .539% |
| 2<r<4 | .0882% | .798% |

Total nonaffine kinetic remainder is1.627003116568822e-7. Calling all of it
thermal would already bias the initial pressure. The prior constant-specific-u
closure would also impose P proportional to reconstructed rho inside a cell,
despite the physical IC's constant pressure. Both issues precede evolution.

Seven-point degree-five triangle quadrature and one subdivision differ by up
to6.86e-4 in E/U_physical and6.52e-4 in local relative mass. These are only two
quadrature levels, not an independent converged integration certificate. The
positive-density fit residual measures a different error. A higher-order
spatial projection needs a separately controlled integration reference.

## Independent kinetic-tensor evolution was proposed and withdrawn

Astra initially proposed evolving Pi=int rho vv^T, defining U=E-tr(Pi)/2,
and a zero-third-central-velocity-moment closure without relaxation. Driver
checked the proposed flux/source before implementing evolution. Exactly:

    f_E - .5 tr(f_Pi) = (U/V) (v_n-w_n)
    S_E - .5 tr(S_Pi) = -P div(v).

The same linear combination survives an augmented HLL flux. Native C++
`check_tensor_thermal_identity_v1.cpp` verifies this for100 anisotropic cases,
with residuals <=4.87e-16 (`tensor_thermal_identity_v1.json`). This is an algebra
audit, not a shock simulation. It exposes that unresolved kinetic energy can
accumulate in a separate reservoir without recovering Euler shock heating.
Re=infinity does not prohibit Euler shock entropy production.

Astra agreed and **withdrew the independent Pi transport decision**. It was
never installed in the hydro integrator. No arbitrary Pi→U damping, energy
projection/deletion or thermal repair is adopted.

## Current next implementation

`ASTRA_EULER_SPATIAL_PLAN.md` supersedes the affine-only evolving model. Reuse
the native CUDA geometry, canonical face fluxes, MPI records/restart and
weight–mass contract in a conservative moving-cell P2 Euler representation.
Carry density, both momenta and total energy with the same spatial order;
internal rotation/strain kinetic energy and angular momentum are derived from
the reconstructed field. Pi/R remain diagnostics. The initial conservative
projection and cold pressure/positivity errors must be quantified before
evolution. Preserve all old outputs and keep the target t6–10 and longest
GIZMO evolution unchanged.
