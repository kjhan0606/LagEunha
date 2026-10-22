# Actual finite P3 material-map RHS and internal rotation

2026-10-21. Native C++/OpenMP/MPI/CUDA. This is an instantaneous RHS at a
deformed control state, not an accepted time step, remap or Kepler evolution.

## State and implementation

The preserved accepted P3/h2 acceleration in
`material_p3_hybrid_v4/mpi1/level1_Kepler_state.bin` defines the finite control
chi=X+h*v0+h^2*a0/2, v=v0+h*a0, h=2e-4. All current coordinates and velocities
are passed to the finite solve; the initialization function does not reset
them. Gravity is evaluated at current chi, not at reference X. The same fixed
reference mass/trace matrices and physical source domain are retained.

Canonical `Exam/gfs_material_map12.h` adds a distinct12-function map API. It
does not feed12 entries into the old5-function API. The stress-work helper
accepts physical symmetric stress, density, gravity and current map. The
reference test alone sets P=P0*J^-gamma as an adiabatic GEOMETRY CONTROL.
This is not a production thermal evolution law or entropy reset.

Actual mass-matrix moments give cell mass, centroid, momentum, K, centered
polar inertia I and intrinsic angular momentum Lint. Krot=Lint^2/(2I) and
Kint-Krot are diagnostics of the actual resolved velocity field. No separate
spin reservoir, relaxation, physical viscosity or heat correction is added.
The intrinsic angular rate includes the moving-centroid term
dLint/dt=dL/dt-c cross dp/dt. Rotation-energy rate includes dI/dt.

Six controls are evaluated: initial Kepler, curved Kepler, uniform velocity
boost, gravity off, linear gravity, and rigid rotation on the same curved
geometry. Current Piola traction and mapped-boundary GCL are evaluated.
The old initial-pressure traction-difference diagnostic is not used to
represent current physical stress in this finite test.

## Completed native bundle

`MATERIAL_FINITE_BUNDLE_20261021.json` records raw native results and the audit.
Each version has24 physical cells,384 internal triangles,2088 broken scalar
DOFs,1865 conforming scalar DOFs, OpenMP2 and MPI1/4 on CPU and CUDA. CUDA MPI4
shares one RTX5000Ada; this is not a four-GPU scaling result. All8 executions
are terminal with exit0 and AUDIT_COMPLETE.

| Quantity, maximum over all modes/controls | v1 | v2 |
|---|---:|---:|
| quadrature | 32x20 | 64x28 |
| positive volume points | 291840 | 817152 |
| numerical cell energy-rate residual / U | 1.19e-11/time | 2.93e-11/time |
| angular-rate residual / 2K | 1.14e-13 | 9.62e-14 |
| intrinsic-angular-rate residual / 2K | 2.44e-12 | 4.41e-12 |
| centered matrix moments vs positive quadrature | 4.63e-14 | 4.44e-14 |
| mass-moment relative difference | 7.25e-15 | 9.68e-15 |
| matrix K vs quadrature / U | 3.86e-10 | 9.81e-10 |
| GCL volume-rate identity / area | 1.49e-13/time | 2.57e-13/time |
| CPU/CUDA force relative difference | 6.79e-15 | 7.16e-15 |
| rigid-rotation physical thermal rate / U | 1.26e-14/time | 5.53e-15/time |
| rigid-rotation remainder / Kint | 3.09e-14 | 2.04e-14 |

Minimum sampled J is0.99999996416, all sampled pressures positive, invalid
point count zero. Sample positivity is not a proof of global injectivity.
Same-backend MPI1/4 states and forces bitmatch in all6 cases. The geometry
packet bitmatches across all four modes at each quadrature. Uniform boost
leaves forces and accelerations bitwise unchanged.

Quadrature refinement changes physical face acceleration by at most2.31e-12
and traction by4.79e-15 across the controls. Actual curved Kepler force
relative change is8.13e-15; gravity-off force is tiny and its relative change
is2.39e-12. Thermal integrals change by less than3.1e-22 in absolute units.

Warm MPI1 actual-curved force/moment evaluation costs0.04776s CPU vs0.00761s
CUDA including transfers for v1; CUDA kernels alone0.001088s. v2 costs0.13265s
vs0.01943s, kernels0.002857s. First-call/module initialization and static
upload cost about0.22s and are separate. Cell aggregation and the reference
dense hybrid solve remain on CPU. These are stage timings, not whole-solver
speedups; the dense reference Schur is not the full71k production design.

## Actual curved packet and exact ray algebra

`curved_map_control.bin` serializes current x/v, original reference triangles,
stable initial-residual basis data and physical parameters. Its magic differs
from the accepted P3 seed and finite-state output; it is not a restart. No
mass matrix is serialized. Original generator sites are advected by actual
chi(q0); there is no centroid or orbit steering.

`curved_map_algebra_v1` and `v2` independently read the coarse/fine packets.
For this selected P3 plus original radial-velocity-residual space, at fixed
reference radius all current x/v are cubic trigonometric functions. With
u=tan((theta-phase)/2), each physical Laguerre power-plane preimage is a
polynomial of degree<=6 divided by positive (1+u^2)^3. This uses arbitrary
current coefficients, not a prescribed future circular orbit.

Both native algebra runs completed:6144 positive interior queries and282624
power-plane comparisons, including nonzero weight differences. Position and
plane errors normalized by triangle size are below2.89e-14; radial/angular
derivatives and J agree to3.7e-15. Actual edge midpoint departure from its
chord reaches2.1043e-6, or5.079e-4 of triangle length. Thus the packet contains
genuine curved geometry. No roots, overlaps, weight solve or cold
re-expression are established by this algebra check.

## Remaining dependency and final scope

Next is actual curved source-to-Laguerre overlap with controlled geometric
uncertainty and physical moment transfer, then admissible cold re-expression.
Composed density-wave/sliding/shock checks and scalable native runtime are
still required before historical t6/8/10 and the longest GIZMO duration.
No new Kepler time evolution has been run from this finite RHS control.
