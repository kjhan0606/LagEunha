# Approved smooth radial-resolution comparison

2026-10-21. User explicitly approved the driver's proposal. One bundled
review by the existing Astra planning agent selected this design. This
native IC experiment can proceed independently of the incomplete curved
material/remap error-ledger research. That research and its failed/partial
certificates remain preserved; their success is not presumed here.

The preceding opinion-only reply added no calculation evidence. There is
no external blocker. The next action is execution of this approved pair.

## Identity and implementation

Authoritative source repository:
`/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`,
HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625` plus existing driver changes.
Work root is the directory containing this document. The supplied backup
working directory is not authoritative. No existing evolution code is
discarded, reset or overwritten.

Both paths use the SAME frozen native C/CUDA/OpenMP/MPI implementation:
`raw_vertex_dedup/v4/eunha2_raw_vertex_dedup.exe`, SHA256
`f6f51995a76861b8438ee6649185b860ceb671c3209097005c7fdf6e13cb250f`.
Original RK4, old velocity limiter, LSF gradients and raw complete geometry
are retained. There is no newly implemented intrinsic-spin evolution in
this particular control. MPI1+OMP2 is the strict failure-witness evolution
path; a separate graded MPI4 one-step control tests the changed layout.
The existing failure-snapshot implementation supports MPI1 only.

Python constructs ICs, launches native executables, and analyzes their
outputs. It does not integrate the fluid equations.

## Fixed initial mesh law

Let h0=24/256=.09375 and rc=.25. The target particle surface number density
is n(r)=h0^-2 [1+63 rc^2/(r^2+rc^2) S(r)]. S=1 for r<=2, S=0 for r>=4,
and S=1-10s^3+15s^4-6s^5, s=(r-2)/2 in between. Thus the target central
length is h0/8, the intermediate length is approximately proportional to r,
and the outer density is unchanged. These are target lengths, not claims
about every actual cell.

The matched uniform control has n=h0^-2. Both use scrambled Sobol rejection
sampling with fixed seed17012026, followed by16 deterministic initial
periodic weighted-CVT steps. In2D the CVT integration density is n^2.
Neither path uses concentric rings. The CVT monitor is distinct from gas
density and from Laguerre weights. CVT is used only before the IC is saved;
centroid steering during native evolution is0.

Physical conditions remain the original smooth density profile, P1e-6,
stored gamma(float1.666667), GM1, eps.047, Rin2/Rout8, no self-gravity,
Re=infinity, W0, individually time-independent masses. Each initial mass
is an integral of the SAME continuous density over its actual periodic
Voronoi polygon. U=P V/(gamma-1). Velocity is the original pointwise
Plummer circular velocity at the generator, as in the preserved coarse/
abrupt pair. K/L differences due to this discretization are reported;
thermal energy is not changed to match totals.

## Preserved historical comparison

The abrupt r<.5 split8x8 case had71080 cells versus65536 and reduced first
positive core Udot about973x. It failed strict RK stages near t.102 (old
limiter) and t.119 (vector limiter). Those data are not rerun or relabeled
as successful. The failure cells were deep in the core, so any improvement
from the new law will not establish that an abrupt interface was the sole
cause. See T8_REVIEW_AND_NEXT_BATCH.md and core_refinement/results_native_v4.json.

## Execution and gates

1. Audit actual area, density integration convergence, total mass/heat,
   initial K/L, radial cell size, adjacent size ratios, inertia eigenvalue
   ratios, polygon widths and neighbor counts. Preserve IDs/initial radii
   and input/source/binary SHA values.
2. Native one-step geometry: positive stages, independent cell volumes,
   surface closure/shared-face reciprocity, absolute force/work residuals,
   and no clipped halo search. Graded MPI1/4 must agree before its longer
   evolution. Keep tiny FP faces visible in the audit.
3. Evolve from saved positive checkpoints through t.05/.1/.2 with identical
   native settings. Strict nonpositive U, any energy-floor injection,
   nonfinite state, geometry/halo failure stops that path and preserves its
   normal step-start and invalid stage. No automatic dt retry or heat fix.
4. If the graded path remains physically accurate, perform ONE effective
   half-step control and then continue its preserved checkpoint through
   .5/1/2/4/6/8/10. Do not repeat completed intervals. Actual accepted dt/RHS
   counts must differ before calling the control a time-resolution test.

Fixed engineering promotion checks in every initial-radius bin below8:
mass-weighted entropy-ratio L1<=.1; .1<=U/Uad<=10 in every cell;
mass-weighted orbital position and velocity errors<=.01. These are early
continuation checks, not published GIZMO thresholds or final accuracy
criteria. A positive but severely depleted/heated solution does not merit
an expensive t10 continuation. A failure of the uniform control does not
automatically reject an independently valid graded path.

Initial radial bins: <.047, .047-.1, .1-.2, .2-.5, .5-1,1-2,2-4,4-8,
plus outer background. Report signed/absolute entropy change, U/Uad,
orbital errors and radial drift, L transport, density-field accuracy,
mass/W invariance, stage/floor evidence and actual times. The frozen face
packet lacks reconstructed face velocities: generator normal jumps are
separately labeled, not passed off as reconstructed jumps.

## Ownership, storage and completion

Preparation is on user-authorized lageunha, not syntax. All executing
sources are frozen before launch. Native Slurm jobs request explicit
memory and `--gres=gpu:H200:1`; independent IC paths may run concurrently.
Two ICs occupy approximately90MiB. Only named checkpoints, first-force
audit packets and failure witnesses are retained; no per-step full dumps.
Estimate remaining wall time from observed accepted dt and step cost.

Passing t.2 is not the goal. The goal remains physically accurate native
CUDA/MPI Kepler across t6/8/10 and eventually the longest GIZMO-paper
duration, including any required physical closure for weights/rotation.
No new long-duration stability has yet been established by this experiment.
