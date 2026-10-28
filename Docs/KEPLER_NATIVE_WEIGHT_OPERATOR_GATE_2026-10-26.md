# Native GPU Laguerre weight operator: implementation gate (2026-10-26)

## Evidence that fixes the target

- C010's first disk ring (`2 <= initial r < 2.25`) reaches mass-weighted
  absolute density changes 0.10488, 0.12553, and 0.14692 at native times 6,
  6.25, and 6.5. Exact-orbit Voronoi reconstructions predict these to within
  3.82e-5, 5.33e-5, and 7.68e-5. This identifies the dominant early disk
  signal as fixed-mass Voronoi geometry under Kepler shear.
- At actual t=6.25 generator positions, 21,436 disk cells can reach their
  initial volumes to maximum relative error 4.85e-12 with power weights.
  Four inner-boundary faces then have `|delta W|/d^2=1.219`, exceeding the
  current native 0.8 guard; the four affected disk generators also fail the
  current own-cell containment ceiling. No static weight cap or simple taper
  tested so far removes this without creating large local density error.
- The signed solution has 83.5% negative disk weights in the zero-exterior
  gauge. A common +0.0036761 shift to **all** sites makes W nonnegative and
  leaves cell areas unchanged, but does not cure those four face violations.

## Operator contract

For a power face between sites `i,j`, with generator positions `x`, velocities
`v`, weights `W`, face centroid `c`, and `d=|x_j-x_i|`, the normal geometric
face velocity is

`u_f,n = [x_j·v_j - x_i·v_i - (v_j-v_i)·c + (Wdot_i-Wdot_j)/2]/d`.

The face-length sum for each cell must match its actual power-polygon area
derivative. An independent t=6.25 calculation confirms this identity to
relative L2 error `2.50e-7`. The static volume Jacobian has each edge
conductance `g_ij=L_ij/(2d_ij)` and gives the stage equation
`J Wdot = target_Vdot - geometric_Vdot_at_fixed_W`.

For the Kepler equilibrium gate, `target_Vdot=0` for active disk cells and
generators follow the existing unsteered motion. A 32-iteration CG reference
reduces the first-ring area change in a `delta t=0.001` Euler geometry check
from median `8.27e-5` to `2.99e-7`, with 95th percentile from `2.72e-4` to
`1.32e-6`. This is one geometry step; do not infer hydro stability from it.

## Native code attachment points

1. The C010 GPU tessellation leaves the directed FaceCSR **on device** before
   the force kernel (`Exam/exam_gpu_extract.c`, `getAccVoro2DBlend_GPU`). A
   host solver would force a face download and repeat host packing. Assemble
   each cell's Laplacian diagonal, fixed-W volume rate, and 16/32-step Krylov
   matvec from the existing device FaceCSR. Download only Wdot per real cell
   if the host RK4 stage update requires it.
2. `Exam/exam_gpu.cu::dev_get2dUpqradRk4` already has the exact Wdot term
   under `w2lin`. Its `laguerre_rot` centroid correction supplies the
   remaining geometric term. The work operator must select the same geometric
   face speed (`geom_fv`) to pass a discrete volume/work audit. Existing
   `GFS_W2_LINEAR` drives a different hydro-kappa target and cannot simply be
   switched on for this solve.
3. The new mode must keep the W=0 Voronoi/fixed-mass reduction, bypass the
   `kappa=0` weight reset only when activated, carry one common nonnegative
   gauge across real/padding/ghost sites, and avoid absolute hydro-derived
   weight clamps. Current `SEDOV_LAGVOL` and CPU-only `GFS_DYN_W2` solve
   different targets.
4. Resolve the four boundary faces explicitly. First test whether the native
   tessellator and conservative pressure/work update remain valid with the
   generator outside its own power cell; positive polygon area alone is
   insufficient. If not, solve a constrained volume problem and use the
   user's approved weight-linked mass exchange with equal-and-opposite face
   fluxes. Do not hide the four cells with a one-sided clamp.

## Release gates before long evolution

- Native H200 static weighted tessellation from the audited t=6.25 checkpoint:
  same site IDs, all 71,080 sites, positive raw areas, candidate completeness,
  disk areas agreeing with the independent regular-triangulation reference,
  and an explicit report for the four boundary faces. Weighted neighbors
  touching the disk span at most three linked cells at this checkpoint; this
  must be checked in the native GPU tessellator.
- One RK4 restart step: identical geometric face speed in volume and work,
  pairwise mass/momentum/energy accounting, no NaN/negative internal energy,
  measured GPU stage time versus C010. Check W=0 bitwise or tolerance-based
  reduction separately. Record whether mass changes, and by which face flux.
- Continue the **same checkpoint chain** through t=10 and past the old t=6--10
  failure window. Then test the paper-matched Hopkins initial condition and
  density field through the published long endpoint; C010's tapered 1/r disk
  is a separate early diagnostic and cannot establish a GIZMO comparison.

Cell spin and angular-momentum exchange remain in the full objective. The
volume experiment shows spin alone cannot remove the dominant fixed-mass
Voronoi density noise; assess spin after the weighted face and mass/energy
operator is conservative, with Kepler physical viscosity kept at Re=infinity.
