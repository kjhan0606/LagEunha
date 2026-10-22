# Shared H1 kinematics: global conservation holds; physical cell ledgers differ

2026-10-21. Continues `ASTRA_MATERIAL_ALE_PLAN.md`. Previous goal turn was
progress (local material RHS, finite affine remap and CUDA helper); this
turn implements the previously absent global continuous kinematics and
consistent mass solve, including native MPI1/4 reference checks. This is
not a full ALE update or a long Kepler stability result.

## Actual implementation

Source remains the confirmed LagEunha worktree at HEAD8e25dfe plus working
changes. New shared native C/CUDA helper: `Exam/gfs_material_triangle.h`.
Reference executable sources: `test_material_shared_v1.cpp`; launcher
`run_material_shared_v1.sh`. They run on authorized lageunha CPU/Intel MPI,
OMP_NUM_THREADS=1, fork bootstrap. No Python hydro evolution.

Use the original24 actual polygons. Exact coordinate keys identify73 shared
P1 vertices/centers and96 positive centroid-fan triangles. The patch has5
connected components (20 adjacent core cells plus4 isolated outer controls).
Each connected component has two continuous input-velocity residual functions

    R_k(X)=v_initial,k(X)-I_P1(v_initial,k)(X).

Their coefficients are global to the component and ordinary kinematic
unknowns, with explicit mass-weighted normalization. The original Kepler
velocity is represented exactly as its shared P1 interpolant plus R. No
initial gravity/derivative-function sweep or orbit reset is used. The scalar
space has83 DOFs for Kepler,73 for the affine shear control. Position and
velocity use the same H1 space; mass is assembled consistently across cells.
Two reference enrichments are kept separate on disconnected components to
avoid artificial coupling across the missing portions of this test patch.

Assemble the same symmetric stress/work quadrature, actual gravity source
and matching-state full Euler boundary traction used by the local test. The
global consistent M solves both components of acceleration. Current stress
is smooth-only sigma=-P I; this is not a shock-stress implementation.

## Shared field, global ledgers and MPI results

Actual same original IC and b=(.013,-.007):

| Quantity | Measured maximum |
|---|---:|
| initial velocity error | 9.88e-16 |
| shared velocity jump | 4.45e-16 |
| shared acceleration jump | 5.69e-14 |
| global mass-matrix solve relative residual | 8.48e-16 |
| global energy-rate residual / sum U | 1.06e-12/time |
| global angular-rate residual / sum 2K | 5.31e-17 |

`material_shared_v1/mpi1` and `mpi4` use the same equations, ownership by
cell modulo MPI size, and explicit Allreduce assembly. Native state arrays
are exported for comparison. See `MATERIAL_SHARED_BUNDLE_20261021.json` for
the actual1/4 differences; normalization of the global residual functions
can change roundoff coefficient values, so bitwise equivalence is not claimed.

## The local physical Euler ledger is not inherited

For each physical cell retain M_i and the true pressure/gravity weak force
f_i. Compute its actual dK_i=v^T M_i a, physical Udot_i, boundary traction
work W_i and same-point gravity power G_i. Also compute the actual angular
derivative, boundary traction torque and gravity torque.

The global equation M a=sum f_i does not imply M_i a=f_i for every cell.
Measured residuals are:

| Case | max cell energy-rate residual / U | max cell angular residual / 2K |
|---|---:|---:|
| exact affine shear control | 3.85e-11/time | 5.84e-17 |
| original Kepler | 15.7337/time | 9.283e-5 |
| original Kepler plus uniform b | 15.9180/time | 9.401e-5 |

These cell residuals cancel globally. The pressure derivative itself remains
at the small initial divergence level: the energy-rate residual above is
**not** a measured heating rate, Pdot, or thermal repair. It is a mismatch
between the globally projected acceleration's cell kinetic derivative and
the physical Euler cell traction/source ledger. Whether the FE method's
implicit intercell transfers can be given an acceptable actual numerical
traction interpretation must be evaluated; this is not by itself proof
that all compatible H1 material methods are impossible or unstable.

## Native source controls identify the discrepancy

`material_shared_v2/mpi1` and `mpi4` keep the exact same input velocity,
geometry, rho, P, H1 basis and mass assembly. Only the diagnostic body force
is changed; these are source controls, not changed Kepler production tests.

| Source control | max cell energy residual / U | max cell angular residual / 2K |
|---|---:|---:|
| gravity off | 6.54e-12/time | 4.23e-15 |
| linear representable g=(-X,-Y) | 6.54e-12/time | 7.43e-13 |
| actual softened Kepler g | 15.7337/time | 9.283e-5 |

The reference kinematic space contains constants/coordinates but not the
arbitrary softened gravitational acceleration. This localizes the effect to
consistent-mass source projection/intercell coupling, rather than to pressure
work or a failure of shared kinematic traces. Adding exact g and its later
derivatives as more benchmark functions is not selected.

Export the actual cell load defect r_i=M_i a-f_i. Its virtual work and torque
reproduce the energy and angular residuals:

    Edefect_i = v^T r_i,
    Ldefect_i = x^T r_iy-y^T r_ix.

For actual Kepler the difference between these identities and the independent
physical ledgers is <=9.95e-11/time when normalized by U for energy and
1.66e-15 when normalized by2K for angular. The large15.73 value is reproduced,
not removed. No added spin, torque correction, heat compensation or source
redistribution is used to conceal it.

## Evaluation dependency and remaining scope

Before promoting this to actual time integration, Astra must evaluate the
local transfer structure: e.g. whether a genuinely derived hybrid/multiplier
interface traction can satisfy the required work, angular, consistency and
entropy contracts, or whether the kinematic/momentum structure must change.
Calling a postprocessed residual an Euler Riemann flux is not sufficient.

Curved Kepler remap, conservative finite re-expression, fixed-target density
wave/sliding tests, shock stress, actual native runtime and t6-10 remain
unverified. Existing local RHS and affine remap controls remain valid in
their measured scope. This shared-field test is a new constraint on the
next implementation, not grounds for starting a long run or claiming full
material-ALE completion.
