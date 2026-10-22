# Compatible material map: local smooth RHS and finite affine remap evidence

2026-10-21. Progress toward the unchanged native Kepler t6-10 and longest
GIZMO objective. This is a **partial first implementation bundle**, not a
completed ALE method or a stability result. No new Kepler timestep or evolved
snapshot was produced. Astra's structural decision is recorded in
`ASTRA_MATERIAL_ALE_PLAN.md`; the previous nodal pair-flux candidates remain
excluded from evolution.

## Source, execution and actual scope

Authoritative source remains
`/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`,
origin `git@github.com:kjhan0606/LagEunha.git`, HEAD `8e25dfe0` plus retained
driver working changes. New native C/CUDA header: `Exam/gfs_material_map.h`.
It supplies actual map/J/velocity-gradient evaluation, symmetric-stress
force/work, mass-matrix kinetic energy, and intrinsic angular/rotational
moments of the represented field. Existing runtime files were not replaced.

All native calculations ran on authorized lageunha, CPU C++ for reference
integration/remap, and its RTX5000Ada for the shared CUDA helper. Python only
reads results; it does not evolve fluid states. Sources, binaries and SHA
checks are frozen separately in:

- `material_rhs_v1`: 24 actual cells, rigid/shear/compression/Kepler/perturbed
  controls, reference quadrature32x20.
- `material_rhs_v2`: finer64x28 Kepler-only comparison with prior native
  initial integrals.
- `material_rhs_v3`: same fine field plus actual M acceleration solve and
  intrinsic rotational moment diagnostics.
- `material_remap_v3`, `material_remap_v4`: finite shear/weight/common-refinement
  controls, with v4 adding derivative convergence and full W/m arrays.
- `material_cuda_v1`: native GPU helper agreement/work/torque controls.

Every native handle is terminal. `MATERIAL_ALE_BUNDLE_20261021.json` gives the
machine-readable summary, hashes, measured limits and incomplete gates.

## 1. Physical internal energy and actual kinematic energy

In material coordinates, use x_h=sum N_a x_a, v_h=sum N_a v_a and
dm=rho0(X)dX. M_ab=integral N_a N_b dm is fixed in the Lagrange stage.
K=.5 v^T M v; physical e0=P0/((gamma-1)rho0) is initialized independently
from the same input field. U is integral e dm, never E-Kaffine.

Each local reference space contains constants and coordinates. For the
actual Kepler IC it additionally contains the two input velocity functions,
after their constant/affine parts are removed stably. The same finite five
scalar functions represent position and velocity. A mass-weighted basis
construction uses the actual Gram, no diagonal ridge or state clipping.
Coefficients are ordinary unknowns; the RHS solves M vdot=force. No exact
Kepler force or residual is subtracted, and no velocity is reset.

Pressure stress sigma=-P I and gravity use the same quadrature and actual
velocity gradient in momentum and physical internal-energy work. Matching
initial boundary traces sample the existing complete lab Euler fan at the
material normal speed; traction and its work are measured, not inferred
from only a pressure average.

Important limitation: this is currently a **cell-local reference RHS**.
Initial traces coincide because the input field is the same, but a shared
global continuous kinematic space/mass assembly is not implemented. Local
accelerations cannot simply be evolved as independent conforming cells.
No new reference-coordinate thermodynamic runtime/checkpoint format exists.

## 2. Actual original24-cell Kepler and perturbation results

Same rho,v,P, eps=.047, gravity, gamma and b=(.013,-.007). No density cusp
smoothing or changed physical IC. Fine64x28 measurements:

| Quantity | Max across both Kepler cases |
|---|---:|
| initial velocity reconstruction absolute error | 8.74e-14 |
| coordinate reproduction absolute error | 1.39e-17 |
| mass-matrix solve relative residual | 2.13e-16 |
| force versus M*acceleration kinetic-power residual | 1.02e-20 |
| energy/stress/gravity work relative residual | 3.65e-15 |
| angular rate residual / sum(M rho v^2) | 2.75e-15 |
| initial abs(div v) | 1.07e-14/time |
| initial abs(Pdot)/P0 | 1.78e-14/time |
| mass difference from prior native integrated IC, relative | 3.40e-13 |
| E difference from prior native integrated IC / physical U | 1.16e-11 |

Changing quadrature32x20 to64x28 changes mass relatively by at most3.09e-15,
and E by at most5.18e-12 of physical U. No nonaffine kinetic energy is assigned
to heat. The earlier affine-only10% initial thermal bias is absent in these
local initial integrals because the actual input velocity field is retained.

Actual material-map volume rate agrees with independent boundary integration
to2.23e-13/area. A centered finite map x+=h v, h=1e-5, gives volume-rate
agreement1.65e-11/area with positive J. This is a geometry difference check,
not evolution to h or proof of later map validity.

The original rigid and simple-shear controls have essentially zero physical
Pdot; compression has Pdot/P0=.500000095/time, matching -gamma div v for
the actual float-stored gamma. v1's relative power metric is ill-conditioned
when both physical work and net boundary work are zero; its absolute errors
are <=7.18e-22. v3 uses absolute virtual-work terms in the normalization.

### Actual internal rotation diagnostic

The mass matrix gives L_int=integral (x-c) cross rho(v-vmean), polar inertia I,
and the rigid rotational projection Krot=L_int^2/(2I). Remaining represented
kinetic energy is diagnostic, not a separate collisionless state or heat.
All measured remainder values are positive in the fine Kepler sample.

Example ID68527, unperturbed:

- L_int=3.04182438e-10;
- internal kinetic energy1.41389542e-8;
- rigid rotational component1.41249679e-8;
- remaining represented strain/curvature kinetic energy1.39863445e-11;
- physical thermal energy2.05993554e-10.

Its actual intrinsic angular RHS is also recorded. This implements the local
torque/energy algebra; it is **not yet a demonstrated evolving cell spin**.
No spin relaxation or finite physical Re was added.

## 3. Finite material deformation and mass-coupled Laguerre remap

Use those same24 initial polygon geometries, with a clearly manufactured
control field rho=1.3, P=1e-6, v=(1.2*y,0), coordinates about(12,12).
The finite material map A=[[1,.125],[0,1]] is the exact simple-shear flow map.
Its J=1, v, physical internal energy, total energy and angular-momentum integrals are
unchanged. Generator positions follow A q; no orbit correction is applied.
The closed test domain is the union of the original24 source polygons, not
an assertion of the full periodic Kepler domain's boundary conditions.

Intersect actual deformed source polygons with target **power support planes**
to form common refinement. Retain the same affine physical field on each
fragment, and integrate m,p,E,L,U,K with positive exact polynomial triangle
moments in long double. This preserves angular information beyond mean p.

Finite target weights solve actual remapped mass, not a rate approximation:

    F_i=m_i(W)+mu_i W_i-m_i0=0,
    J_F=diag(mu)+weighted power-face Laplacian,
    edge coefficient = integral_face rho ds/(2d).

The face measure uses the union of source-domain intervals, avoiding double
counting coincident internal source faces. mu_i=m_i0/V_i0>0. Dense Cholesky
and Newton solve this small reference control; this is not a production
CUDA/MPI finite-weight algorithm.

| Control | Newton | max constraint residual/m0 | max abs W |
|---|---:|---:|---:|
| identity W0 | 0 | 0 | 0 |
| shear to actual Laguerre target | 2 | 4.97e-14 | 2.102e-7 |
| inverse deformation and return remap | 2 | 1.28e-13 | 5.29e-18 |

Both finite solves require no backtracking. Across the forward/return:
mass error <=2.04e-20, px error0, total E error <=9.54e-12 of physical thermal
energy, angular error <=1.09e-18, thermal relative error <=2.28e-19. No kinetic
projection loss was put into U. Sum(mu W) is about1e-22.

The nondegenerate finite target's analytic Jacobian matches independent
central differences at2.08e-13/mu at the first interval. At cocircular W0
and the returned geometry, new diagonal faces activate: the measured error
8.89e-6 ->4.45e-6 ->2.22e-6 halves with interval. Do not label these latter
finite-interval checks machine-precision derivatives at a topology change.

Uncompressed fragments grow24 ->92 ->180. This is allowed only as an oracle
to separate geometry transfer from representation compression. Merely
retaining all fragments is not a scalable long-time method. Replacing each
target field by only mean velocity would lose about424 times the control's
thermal energy in K and lose angular information; this loss is measured
hypothetically and **not applied**. Affine fields remain exact here; curved
Kepler compression has not been tested.

### Corrected diagnostic implementation issues

v1/v2 rejected a -0 fan area. They built overlap clipping normals from rounded
target vertices, creating fake slivers at coincident interfaces. v3 clips
directly with the original power supports; identity produces exactly24 source
pieces and all positive-area fragments are retained. Exact zero-area
face/vertex-only intersections contribute zero. No finite area/length cutoff,
negative-area absolute value or mass fallback was introduced.

The raw v3/v4 JSON has a serialization typo `.125`; it is preserved verbatim.
`material_remap_v4/analysis.json` explicitly normalizes that one literal to
`0.125` in its derived parsed record. The source producer is corrected for
future runs. No physical numeric value is changed or simulation rerun merely
to fix the serialized spelling.

## 4. Native CUDA core check

`material_cuda_v1/contract.json`: 65536 manufactured material-map points on
NVIDIA RTX5000Ada, invalid0, minimum J=.9857569.

- Max CPU/GPU difference2.23e-16.
- Stress/kinetic/thermal/gravity work identity2.13e-16.
- Symmetric internal-stress torque identity1.16e-16.
- Single CPU thread6.21ms, GPU kernel .284ms.

The timing excludes allocation, transfers, basis/mass assembly, remap, MPI
and the real solver. It measures this helper only; no full Kepler acceleration
factor or performance claim is made. The old native CUDA/MPI runtime remains
intact and has not silently become this standalone reference program.

## Required next work before any actual evolution

1. Shared continuous kinematics/global mass assembly; local matching initial
   traces alone do not keep later material faces joined.
2. Actual nonaffine curved material geometry, remap of the original
   nonuniform Kepler field, and cold-conservative finite re-expression.
3. Fixed-target remap density-wave response of the composed update, plus
   the historical cs/h sliding-mode screen; no favorable material-frame-only
   substitute.
4. Multidimensional symmetric compression stress derived from complete Euler
   Riemann data, force/work heating, positivity and Sod/Sedov. Current stress
   is smooth-only tau_num=0; matching-state fan checks are not shock tests.
5. Actual native CUDA/MPI state/halo/migration/restart/timestep integration,
   then the historical t6-10 interval and longest matched-GIZMO duration.

The first three remain parts of Astra's unfinished initial bundle. This
report does not authorize a long run based only on the successful local
operator and affine-remap controls.
