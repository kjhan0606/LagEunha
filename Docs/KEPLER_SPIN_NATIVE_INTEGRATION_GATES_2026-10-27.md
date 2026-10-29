# Native cell-spin integration: required discrete contracts

The user's target includes spin, evolving weights when useful, no generator
steering, fixed individual mass at `W=0`, and inviscid Kepler flow
(`Re=infinity`). The existing native CUDA/MPI/RK4 production path currently
has no evolved spin. The two-cell ledgers and affine reconstruction are
isolated prototypes, not native physics. The following contracts identify
what must be implemented and tested before a spinning Kepler result can be
claimed.

There is already a distinct geometric face-speed correction in the source.
For `W=0`, the Voronoi face-centroid correction is applied automatically on
both CPU and GPU paths. `GFS_LAGUERRE_ROTATION=1` extends the same geometric
correction to nonzero-weight Laguerre faces; the active W=0 long-run scripts
in q1 and 420² leave that extension disabled, but still take the automatic
Voronoi branch.
This corrects the face-normal ALE/Riemann speed and is not an evolved cell
spin: it adds no intrinsic angular-momentum state, torque transport, or spin
energy reservoir. Keep these two mechanisms separately identified in
diagnostics and claims.

## State and initialization

For each cell define mass `m`, momentum `p`, **mass centroid** `C`, total
energy `E`, intrinsic angular momentum `S`, and a central mass-moment
tensor `Q` with `I=trace(Q)`. If an affine velocity representation is
used, distinguish its non-rigid kinetic part `K_nonrigid` from thermal
energy `U`:

```
E = |p|^2/(2m) + S^2/(2I) + K_nonrigid + U + (defined potential term).
```

The initial total-energy convention, potential accounting and equation of
state must be explicit. Adding reconstructed spin by subtracting it from the
current cold-disk `U` is not generally viable: the
[saved-state audit](KEPLER_SPIN_SHEAR_LOCAL_MODEL_2026-10-27.md) finds many
cells would become nonpositive under that hypothetical conversion. A new
consistent energy reference is a possible alternative, but needs its own
initial-condition and conservation tests. The generator position must not
be silently substituted for `C`; first mass moments and their volume/face
transport need a closure. A uniform polygon density identifies `C` with its
area centroid, which may conflict with zero face mass flux on a moving mesh.

## Pair balance and physical closure

For a physically specified shared traction `F` applied to cell i and `-F`
to j at face point `c`, fixed masses and centroids moving with their bulk
velocities give torque increments

```
dS_i = (c-C_i) x F * dt
dS_j = -(c-C_j) x F * dt.
```

These cancel the pair's orbital angular-momentum change. For nonzero mass
transfer, use a consistent first-moment and angular-momentum flux; the
[pair ledger](gfs_spin_affine_face_exchange_20261027.h) is an algebraic
testable example, but does not derive the flux. Use one common, equal-and-
opposite total-energy face flux and obtain `U` after subtracting the other
energy components; ensure positivity without silently injecting energy.
At `Re=infinity`, pressure traction is normal to the face and no physical
tangential viscosity or friction is introduced. A numerical Riemann solver
can have its own dissipation and must be measured, not called physical
viscosity by assumption.

The existing W0 variational pressure operator has a globally conservative
generator force and geometric work, but its per-face traction decomposition
is not unique. Therefore choose and justify a specific pressure/traction and
energy flux before deriving a torque; do not attach an arbitrary torque to a
gauge-dependent decomposition. If spin alters a face velocity reconstruction,
remove the rigid component already present in the existing gradient before
adding the evolved spin. The [affine operator](gfs_affine_spin_reconstruction_20261027.h)
shows that algebra, but does not define how `K_nonrigid` evolves or how a
Riemann state is obtained.

At `W=0` maintained in time, each mass must remain fixed and generator
motion must not be regularized by steering. If nonzero weights are allowed
to change mass, define the actual weight-linked mass, momentum, energy and
moment fluxes and show they reduce to the approved W0 condition. The tested
global projection of weight work is not by itself a proof of local
moving-face transport. Turning the spin feature off must preserve the native
baseline; enabled `S=0` is a separate reduction that must be specified and
tested because the original gradient may already contain rotation.

## Native integration and acceptance

1. Implement a versioned, opt-in state and restart format, including any
   centroid/first-moment, spin and non-rigid state; preserve existing
   checkpoints and MPI ownership/halo identity. A restart must reproduce a
   continuous run at the same temporal stage.
2. Compute required polygon moments in O(number of faces), without global
   subcell quadrature; stage them consistently with the native CUDA/MPI RK4
   geometry and forces. Apply pair fluxes, torque and energy to the same RK4
   stage and ensure race-free GPU/MPI accumulation and rollback.
3. Independently verify mass, linear/angular momentum, total energy, thermal
   positivity, zero unintended injection, and spin/non-rigid work on analytic
   rigid rotation, affine shear, translating cells, equal-pressure and shock
   controls. Establish the disabled and W0 fixed-mass reductions and actual
   restart equivalence before any Kepler claim.
4. Compare matched W0 and evolving-weight, no-steering, `Re=infinity`
   Kepler runs at least through the historical `t=6–10` failure window,
   with density fields, orbital drift, entropy and conservation. Continue
   only on meaningful gates toward the GIZMO figure's `t=600` benchmark;
   an aligned bounded control or a short snapshot is not a substitute.
5. Measure full-step overhead and scaling on the user's native GPU/CUDA+MPI
   solver against a matched baseline; a single-cell device test or arithmetic
   operation count does not establish AREPO-like throughput. Extend spin from
   a scalar in 2D to vector angular momentum and a 3D inertia tensor before
   claiming a three-dimensional method.

## 2026-10-28: first native torque candidate

`W=0` is a geometric statement: every Laguerre weight vanishes, so the
power cells become Voronoi cells. Fixed mass is a separate solver invariant;
zero weights alone do not imply zero mass flux. Keep checking mass by particle
identity whenever the weight or face-motion path changes.

For the first spin-enabled native candidate, take torque from one shared
inviscid ALE face flux, rather than assigning a moment to the existing
variational generator force. The latter has a non-unique reciprocal
face-force decomposition, so its generator force does not determine a unique
cell spin source. With outward normal `n`, face area `A_f`, face centroid
`x_f`, and one common Riemann pressure `P*`, use

```
F_i = -P* A_f n
tau_i = (x_f - C_i) x F_i
```

and the equal-and-opposite force on the neighbor. For a zero-mass-flux face,
the same face state must use the pressure-work flux `P* u*_n A_f`; including
the reconstructed rigid and non-rigid velocities in `u*` makes its power
consistent with bulk, spin, and affine kinetic energy. If mass flux is enabled,
the momentum, energy, first-moment, and angular-momentum fluxes must all use
the full ALE Euler flux and reduce to the fixed-mass expression when the mass
flux is zero. This candidate adds no physical viscosity or spin relaxation
at `Re=∞`.

This is a physics selection for a separate opt-in branch, not a claim that the
current variational W0 path already implements it. The implementation still
needs a shared-face GPU/MPI/RK4 update, a consistent mass-centroid and `Q`
evolution, extended-energy initialization, and the native reduction and
restart gates above. Compare its disabled and enabled trajectories against
the unchanged baseline before any long-time Kepler run.

## Active-source audit (2026-10-28)

The active source directory is `src/`, at commit `75aaeed`; it also contains
existing uncommitted edits in `Exam/exam.c`, `Exam/exam_gpu.cu`, and
`Exam/exam_gpu.h`. These were preserved during this audit. The separate
polygon-moment candidate under `/tmp/Kepler_spin_moments_20261028` is based on
another repository state, so its patch must not be copied wholesale into
`src/`.

The active `RK4` record in `src/eunha.h` stages position, velocity, internal
energy and entropy, while the RK4 updates are explicitly assembled in
`src/Exam/exam.c`; there is no spin stage yet. `ParticleSoA` in
`src/Exam/exam_gpu.h` has no mass-centroid, central-moment or spin arrays, and
its `PolygonOut` currently stores polygon vertices and connectivity only.
The standalone `RunKepler` driver uses raw `kpout` particle records. The active
q1/q2 experiments instead use `RunKH` and `khout`: its current header stores a
magic value, particle count, record size, and stress flag, then writes complete
stress-particle records. The `khout` reader already has a common-prefix path
for legacy/smaller records. Affine state therefore has to extend this existing
`khout` versioning and remain in the MPI-migrated stress particle; a separate
sidecar would lose identity/state when particles migrate.

## Active-job driver and force-path correction (2026-10-28)

The earlier audit followed the standalone `RunKepler` function and incorrectly
assigned that path to the active Slurm experiments. Direct inspection of jobs
437771 and 438040 shows the active executables use the KH driver with the
Kepler IC/force enabled by `EUNHA_IC=kepler`. Job 437771's executable hash is
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`; its ELF
symbols include `RunKH` and `exam2d_vph_rk4_int_blend`, and omit `RunKepler`.
Its startup log reads `khout.000014.dat`, reports `av_mode=5 use_stress=1`,
and initializes an A10 GPU. The q1 job has the same KH/blend symbols and also
uses the GPU path.

The binary's debug source path resolves to
`/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`,
at commit `472e7fcc`; `Exam/KH/kh.c` selects
`exam2d_vph_rk4_int_blend` for `use_stress && GAS_EVOLMETHOD==1`. That RK4
path contains the central Kepler acceleration and dispatches force evaluation
to `getAccVoro2DBlend_GPU` when CUDA is enabled. The standalone
`RunKepler`/`exam2d_vph_rk4_int` CPU path remains a separate driver and is not
the path used by these active q1/q2 jobs.

Therefore the first integrated affine-spin implementation for the current
experiments must target the stress-particle blend RK4 and its CUDA face-force
path, including device state transfer and MPI migration. The `src/` worktree
at `75aaeed` remains the current development tree with user changes; its
build must be matched to the running-job baseline before interpreting any
comparison.

## 2026-10-28: conservative affine cell state

A torque-only scalar spin update is insufficient: the current `ie` stage gets
the pressure-work residual after bulk kinetic work, so adding rotational work
without changing that residual can double-count or omit energy. The cell state
should carry conservative first, second, and momentum-cross moments together
with total energy. This represents rigid rotation and affine shear in one
velocity field; its antisymmetric part is the cell's intrinsic spin. It adds no
viscosity or spin damping at `Re=infinity`.

For one moving control volume, define

```
M       = integral rho dA              P_b = integral rho u_b dA
H_a     = integral rho x_a dA          Q_ab = integral rho x_a x_b dA
R_ab    = integral rho x_a u_b dA      E   = integral (rho e + rho |u|^2/2) dA
```

Use one common Riemann state on each face, with outward normal `n`, face speed
`w`, and `Fmass = rho* (u* - w) dot n`. The ALE Euler moment equations are

```
dM/dt     = -sum_faces integral Fmass ds
dH_a/dt   = integral rho u_a dA
            -sum_faces integral x_a Fmass ds
dP_b/dt   = -sum_faces integral [rho* u*_b (u* - w) dot n + P* n_b] ds
            + integral rho a_b dA
dQ_ab/dt  = integral rho (x_a u_b + x_b u_a) dA
            -sum_faces integral x_a x_b Fmass ds
dR_ab/dt  = integral [rho u_a u_b + P delta_ab] dA
            -sum_faces integral x_a [rho* u*_b (u* - w) dot n + P* n_b] ds
            + integral x_a rho a_b dA
dE/dt     = -sum_faces integral [E* (u* - w) dot n + P* u* dot n] ds
            + integral rho u dot a dA
```

For the fixed-mass reduction, explicitly set `Fmass=0`; the same equations
still evolve the centroid and internal moments under the cell motion and face
stress. `W=0` determines the Voronoi geometry but does not itself set the mass
flux to zero. For nonzero/evolving Laguerre weights, retain the mass flux and
use it in every moment equation. A face flux must be computed once and applied
with opposite sign to its two cells to conserve global mass, momentum, and
angular momentum.

**Qualification added 2026-10-29:** simply imposing `Fmass=0` is not by
itself a valid Euler ALE flux on a generator-defined Voronoi face. The
continuum flux is `rho (u_n - w_n)`, so for positive interface density it
vanishes only when the face's geometric normal speed `w_n` agrees with the
normal material/contact speed. In a moving Voronoi mesh, its speed is
determined by the generators and can differ from the fluid speed; the
[Q2 diagnostic](KEPLER_Q2_FACE_SLIP_2026-10-29.json) provides a kinematic
example, although it does not itself measure the Riemann contact speed.

For a simple counterexample, take a uniform stationary fluid on either side
of a face moving right with speed `w_n>0`. The Euler ALE mass flux is
`-rho*w_n`; zeroing it while retaining the same moving face is a different
discrete model. Therefore the fixed-mass spin branch needs a derived
non-Euler moment closure or a material-interface representation compatible
with the generator geometry; it cannot claim that the displayed continuum
ALE equations remain exact after zeroing their mass flux. No generator
steering is allowed. The compatibility of that closure with fixed mass,
pressure work, moments, energy and the W0 Voronoi limit is a prerequisite,
not a post-hoc conservation test. A branch that moves mass with evolving
weights also needs the actual geometric face speed and the full consistent
fluxes; weight changes alone do not establish this compatibility.

The running Q1 native log identifies `GFS_VARPRESS_W0`, whose
[CUDA force/work branch](../Exam/exam_gpu.cu) uses a fixed-mass variational
`B^T P` force and `-P Vdot` internal-energy work. It is not the candidate
Euler ALE moment formulation above, and its successful fixed-mass bookkeeping
does not supply the missing moment/traction closure. This distinction also
means a pressure-work-only or shared-normal-traction substitution cannot be
claimed as a spin implementation. Existing W0 runs remain unchanged.

Recover the centroid and central moments by `C=H/M`,
`S=Q-H H^T/M`, and `B=R-H P^T/M`. Under a linear within-cell velocity closure,
`u(x)=P/M + G (x-C)` with `G=B^T S^{-1}`. Then
`E_aff=0.5 tr(G S G^T)`, `E_bulk=|P|^2/(2M)`, and
`E_thermal=E-E_bulk-E_aff`; in 2D the intrinsic angular momentum is
`Lz=B_xy-B_yx`. This gives a direct positivity check for the thermal residual
and avoids an independent spin-energy update inconsistent with pressure work.

With the current midpoint face flux, segment geometry moments are inexpensive
and exact for a constant face state: for endpoints `a,b`, length `l`, midpoint
`m`, and `d=b-a`, use `integral(1)=l`, `integral(x)=l*m`, and
`integral(x x^T)=l*(m m^T+d d^T/12)`. This is O(number of faces), not a cell
interior quadrature.

In the active RK4 path, this state must be staged at all four RK stages and
included in the versioned checkpoint. The existing `ie` pressure-work update
must be replaced by the total-energy residual above when this mode is enabled.
Acceptance requires positive-definite `S`, positive `E_thermal`, restart
equivalence, and conservation of total orbital plus intrinsic angular
momentum and total energy before applying the mode to the Kepler disk.

### Host-side affine-moment prototype (2026-10-28)

An isolated development worktree at
`/tmp/Kepler_affine_spin_20261028` was based on the current active source
commit `75aaeed`, with its existing tracked local changes applied. The user's
`src/` worktree was left untouched. `Exam/kepler_affine2d.h` now contains
double-precision exact polygon area/centroid/central moments, exact straight
face first/second moments, and recovery of the affine velocity gradient,
intrinsic angular momentum, and bulk/affine/thermal energy split. It also has
fixed-mass Lagrangian cell sources and a one-call shared pressure-face update
that applies action/reaction, matching pressure work, and the corresponding
intrinsic moment to both cells. Its initialization helper also computes the
exact gradient of the softened Kepler IC field used by `kp_mkinitial`.

A standalone C99 analytic check compiled with `-Wall -Wextra -Werror` and
passed for a square in both windings, a rigid rotation with `Omega=3`, and a
straight-face moment. It also verified zero shape change for rigid rotation
under uniform pressure, equal/opposite face momentum and energy rates, and
the face torque sign; the softened Kepler velocity gradient agrees with
centered finite differences. The helper is not yet called by the active KH-driver
RK4/GPU path; these checks validate its algebra only. Native RK4 stages,
versioned `khout` restart integration, and GPU/MPI state transfer remain
outstanding.

## 2026-10-29: independent affine moment balance check

The actual `Exam/gfs_moment.h` source operators now have a separate
[analytic test](test_gfs_moment_affine_balance_20261029.c). Unequal point
masses, an anisotropic distribution and a general affine velocity field
independently reconstruct the velocity gradient and thermal energy; direct
point sums verify all first/second moment and momentum-moment volume source
terms, including isotropic pressure. A shared ALE face flux with **nonzero
mass flux** and distinct cell reference origins verifies pair mass, momentum,
energy, first-moment and total angular-momentum balance. A central-gravity
point checks zero torque and consistent gravitational work. On lageunha the
C99 test passed with `gcc -O2 -Wall -Wextra -Werror`.

These are small host algebra checks on the existing moment helpers; they do
not execute the CUDA face solver, RK4 stages, mesh motion, MPI migration,
restart, or a Kepler evolution. In particular, testing a general face mass
flux does not alter the required fixed-mass `W=0` reduction. The production
spin/affine integration and its conservation and performance gates above
remain open.

## 2026-10-29: optional density-fit warm start

The positive exponential-quadratic density reconstruction currently starts
its five-parameter Newton iteration at zero. `gfs_moment_density_fit_warm`
now accepts a finite initial guess in the current polygon's normalized basis,
uses the same moment residual tolerance and line search, and retries the
unchanged zero-start solve if the guess fails. A nonfinite seed also takes
the original path. No state projection, density clipping or thermal-energy
injection is introduced; the original `gfs_moment_density_fit` entry point
continues to use the same zero initialization. This optional API is not
connected to initialization or an evolution path.

In the [host test](test_gfs_moment_density_warm_20261029.c), two nearby
realizable nonuniform states on an asymmetric polygon produced matching
cold/warm fitted parameters within 1e-9, residuals below 5e-12, and **4 vs
2 Newton iterations**; nonfinite seed fallback and rejection of an invalid
moment state passed. A minimal CUDA kernel calling the warm API compiled with
`nvcc`, but was not executed on the GPU. This one favorable synthetic case
does not measure full-step speed or worst-case retry overhead. A moving polygon
changes the normalized basis, so previous parameters are only a guess; actual
stage-by-stage trajectory, positivity, bitwise or tolerance-based baseline
equivalence, device runtime and adverse-seed costs must be measured before
using this option in a solver. It does not make the nonlinear fit itself an
AREPO-speed method.

## 2026-10-29: tentative moment RK4 stage storage

The existing tagged moment payload aliases the legacy particle RK4 scratch;
the legacy integrator must not use both at once. A separate opt-in
[moment RK4 helper](../Exam/gfs_moment_rk4.h) now stores base, tentative stage,
and weighted stage sum outside that legacy scratch. It implements classical
four-stage RK4 without overwriting the accepted state. Each tentative stage
requires positive-definite mass moments, positive thermal residual, and
finite derivatives; an optional fixed-mass contract rejects nonzero mass
RHS. Failure leaves the previous stage scratch unchanged. A completed result
remains tentative: a native implementation must perform an all-cell/all-rank
decision before committing any particle.

The [host check](test_gfs_moment_rk4_20261029.c) reports oscillator error
ratio **15.995** when halving the step, exercises stage rejection and
unchanged input, and checks mass, momentum, energy and angular momentum for a
two-cell conservative exchange. The separate [CUDA check](test_gfs_moment_rk4_gpu_20261029.cu)
compiled for sm_89 and ran successfully on Ada for one cell and ten analytic
oscillator steps. This verifies RK4 algebra and the helper's fail-closed
behavior, not the actual pressure or spin dynamics. Face/geometry stages,
reconstruction, rollback of whole-system side effects, MPI collective
acceptance, KH restart and production integration remain unimplemented.

### Moment-only MPI commit gate

The separate [MPI commit helper](../Exam/gfs_moment_rk4_mpi.h) checks each
tentative cell's RK4 completion, unchanged base, realizability, fixed-mass
constraint when requested, and identical positive `dt` across nonempty ranks.
It commits moment arrays only after an all-rank success vote. Every rank,
including an empty rank, must call the collectives in the same order. The
[two-rank test](test_gfs_moment_rk4_mpi_20261029.c), run with Intel MPI on
lageunha, passed intentional one-rank rejection, different-rank timestep
rejection, successful shared-face conservation, stale-base rejection and an
empty-rank vote. Rejected accepted-state arrays stayed unchanged.

This is still a **moment-array-only** transaction and is not wired into the
native RK4 driver. Its vote does not establish matching physical start times,
unique cross-rank faces, correct halo exchange, or stage-by-stage collective
ordering. Geometry, generator position, pressure/thermal state, face ledgers,
wall events and device buffers require a joint rollback and common acceptance
before it can be used in the real solver. Testing a single shared flux on two
local CPU ranks does not establish CUDA+MPI correctness or a Kepler solution.

The host test also covers an affinely changed polygon, passing the seed
in-place, and a finite but extreme seed that makes the seeded Newton attempt
fail. In that case the warm wrapper recovers the exact cold-start result;
nonfinite seeds similarly fall back. This is a safety check on a few synthetic
states, not a performance qualification for a moving native mesh; a failed
seed incurs the cost of both attempts.

## 2026-10-29: weight-linked mass-exchange building block

An isolated [weight exchange helper](../Exam/gfs_weight_mass_exchange.h)
implements a candidate incremental mass flux associated with the weight-rate
contribution to power-face speed, relative to a separately defined fixed-mass
baseline. With the normal from i to j, generator separation `d`, signed power
weights `W_i`, and one common face state, its contribution is

```
delta_w_n = (Wdot_i - Wdot_j)/(2 d)
F_m       = -rho_face delta_w_n
(F_momentum, F_energy) = F_m (u_face, e_total_face)
```

The helper applies this flux once at one common face point, with opposite
signs and each cell's own fixed moment origin. `Wdot_i=Wdot_j=0` gives zero
exchange, including the approved W=0-for-all-time fixed-mass case. The
[C99 host test](test_gfs_weight_mass_exchange_20261029.c) checks the sign,
swapped-face antisymmetry, mass/momentum/energy, global first/second moments,
angular momentum, Galilean transformation and invalid inputs. It passed with
`gcc -O2 -Wall -Wextra -Werror`; a minimal CUDA caller compiled with `nvcc`.

This is **not a complete Euler ALE flux** or a validated fixed-mass spinning
scheme: the unweighted geometric/material slip discussed above is deliberately
not represented by this incremental term. It also does not supply a physical
face density/velocity/energy state, pressure traction and work, a weight
evolution law, positivity, finite-time geometric remapping, topology-change
handling, race-free GPU/MPI ownership, or RK4 staging. Exact pair conservation
does not prove consistency, convergence, entropy behavior or Kepler stability.
The helper is unused by the native solver and the ongoing runs; those physics
and implementation gates must be closed before activation.

### Independent geometry check

The actual C flux helper was loaded through a [small bridge](gfs_weight_exchange_bridge_20261029.c)
and compared to centered finite differences of independently rebuilt power-cell
areas at fixed generators, with unit density. The geometry was taken from
the saved Q2 t=47.50031653561342 state; all **3,931** inner-disk cells and
their touching faces were checked with a reproducible nonuniform weight-rate
field. For W=0, the relative RMS difference between the flux-predicted area
rate and the finite-difference rate was **2.27e-9** for step 1e-4. Repeating
with small synthetic nonzero **signed** power weights gave **1.82e-9**.
Results for two smaller differencing steps show increasing roundoff, not a
measured convergence order. See [W0 receipt](KEPLER_WEIGHT_EXCHANGE_GEOMETRY_2026-10-29.json),
[weighted receipt](KEPLER_WEIGHT_EXCHANGE_WEIGHTED_GEOMETRY_2026-10-29.json),
and their matching audit scripts.

This validates the infinitesimal weight-speed sign and area coefficient for
these frozen geometries and sample rates; the nonzero signed weights are a
mathematical power-diagram test, not evidence that all native weight paths
support negative values. It does not validate the unresolved material-slip
closure, finite-step conservation through topology changes, physical face
states, pressure work, or a native GPU/MPI evolution.

### Segment-integrated moment correction

The initial exchange helper treated the whole face flux as if concentrated
at its midpoint. Although that preserves global pair sums, it does **not**
give the correct local second spatial moments of a finite segment. The point
interface is now explicitly named `gfs_weight_exchange_point` and takes a
quadrature weight; the new `gfs_weight_exchange_segment` uses two Gauss
points for a constant shared flux, integrating straight-face mass moments
through degree two exactly. A changing reconstructed face state must instead
be evaluated consistently at the quadrature points; a constant-state segment
call does not provide higher-order flux reconstruction.

The updated host test checks the analytic local `xx`, `xy`, and `yy` segment
moments, pair conservation, invalid-segment rejection, and Wdot=0 reduction;
it passes with strict C99 warnings. A minimal CUDA caller compiles but has
not run on a GPU. An independent finite-difference check of power-polygon
area and first/second moments for all 3,931 inner Q2 cells with synthetic
nonzero signed weights agrees with the actual C segment operator to relative
RMS **1.82e-9–2.18e-9** across `M`, `H_x`, `H_y`, `Q_xx`, `Q_xy`, and
`Q_yy` at differencing step 1e-4. See [receipt](KEPLER_WEIGHT_EXCHANGE_MOMENTS_2026-10-29.json)
and [reader](audit_kepler_weight_exchange_moments_20261029.py).
This remains an infinitesimal frozen-geometry test, not proof of finite-time
conservative transport, native negative-weight support, or a working spin
solver.

### Experimental CUDA face accumulation

The isolated [CUDA kernel](../Exam/gfs_weight_exchange_gpu.cuh) evaluates the
same segment-integrated weight-only contribution for a caller-supplied list
of unique interior faces and accumulates all 13 moment derivatives using
double-precision atomics into a disposable cell RHS. A bad face increments an
error counter; because other faces may already have contributed, the caller
must discard the **entire** RHS on any error. The kernel cannot certify that
the face list is unique or complete, and does not handle MPI ownership,
ghosts, wall faces, physical face-state reconstruction, full ALE transport,
weight evolution or positivity. It is not connected to the production CUDA
RK4 path.

The [CUDA test](test_gfs_weight_exchange_gpu_20261029.cu) was compiled for
sm_89 and ran briefly on the lageunha Ada GPU. Six faces concurrently update
four cells; all 13 resulting components agree with the serial CPU reference,
and total mass, momentum, energy and angular momentum cancel. Zero weight
rates give zero RHS and an out-of-range neighbor is reported invalid. The
existing long-running solver remained live. This is a small kernel correctness
test, not an end-to-end solver, MPI test, GPU performance measurement, or
proof that the moment closure suppresses Kepler instability. Atomic ordering
is not deterministic; a native integration needs conservation tolerances and
single-face ownership tests under realistic parallel decomposition.

An isolated [microbenchmark](bench_gfs_weight_exchange_gpu_20261029.cu)
uses a synthetic 235x235 triangular-neighbor grid: **55,225 cells, 164,736
unique faces**. On the same Ada GPU while the native Q1 calculation was
running, ten timed kernels after three warmups had median **0.139 ms**
(range 0.136–0.141 ms). The result's 13 components agreed with serial host
accumulation to relative RMS **1.39e-16**, maximum absolute 2.71e-20;
global mass, momentum, energy and angular-momentum residuals were below
4e-19 in that synthetic setup. The experimental AoS face input occupied
21.1 MB and the RHS 5.74 MB.

This is a short GPU-event measurement with a shared device, not an isolated
benchmark. It excludes construction of unique faces, zeroing, host/device
transfer, EOS/density reconstruction, RK stages, topology changes, MPI,
pressure work and a native full-step baseline. A regular synthetic grid
does not reproduce actual Kepler face contention or nonuniform states. It
therefore cannot establish an AREPO-comparable evolution rate or production
overhead; it only rules out an obvious cost explosion in this one kernel
under the stated setup. The concurrent native process was confirmed live
afterward.

## 2026-10-29: inviscid material-face pressure building block

The [material-pressure flux helper](../Exam/gfs_moment_material_pressure.h)
returns the fixed-mass inviscid face flux `(0, P*n_x, P*n_y, P*u_n)` only
when the caller supplies a unit normal, finite nonnegative pressure, and a
material normal speed matching the **actual geometric** face-normal speed
within an explicit absolute tolerance. A rejected input does not alter the
output. It supplies neither a Riemann state nor a justification for a chosen
tolerance, and cannot turn an arbitrary Voronoi face into a material face.
The normal traction has no added physical tangential viscosity; its off-center
moment can still change intrinsic angular momentum.

The [host test](test_gfs_material_pressure_20261029.c) integrates all four
faces of a uniform-density square with two quadrature points per edge.
Stationary and rigidly rotating affine states show the expected pressure
cancellation, pressure work, moment evolution and zero net torque. A shared
two-cell face verifies equal/opposite momentum and energy and total angular
momentum; a deliberately nonmaterial face is rejected without altering its
flux. It passed C99 `-O2 -Wall -Wextra -Werror`; a minimal CUDA caller also
compiled. These are material-boundary analytic checks only, not an evolved
spin simulation or native CUDA/MPI integration. In particular, Q2's measured
geometric/reconstructed slip means this helper cannot simply be switched on
for the existing W0 Voronoi trajectory. A different consistent closure or a
material-interface representation remains necessary.
