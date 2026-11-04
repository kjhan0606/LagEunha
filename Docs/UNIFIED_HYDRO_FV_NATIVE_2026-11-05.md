# Native dimensional FV--ALE implementation and first evolution results

This is the implementation checkpoint following
`UNIFIED_HYDRO_FV_ALE_PLAN_2026-11-05.md`, not adoption of a qualified new
production solver. The common native C/OpenMP operator now evolves actual
1D intervals and the project's native 2D/3D power cells. Python schedules
executables and assesses saved output; it does not evolve the fluid.
The existing production CUDA/MPI trajectory and accepted Kepler/manuscript
results retain their original state interpretation.

## Implemented scope

`Exam/gfs_fv_native.{h,c}` and `Exam/gfs_fv_geometry.c` implement the shared
five-slot payload `[m,px,py,pz,E]`, actual-cell centroid reconstruction,
fixed scalar/velocity limiters, positive actual-face quadrature, complete
frozen HLLC/HLL ALE fluxes, and the screened weight-rate solve. Birth
`mu=m0/V0^(2/d)`, m0, W0 and the integrated true residual are retained in
each cell. Mass is advanced by its actual flux, not replaced from the
mass--weight relation. Inactive momentum components remain zero.

`Exam/HydroCore/fv_native_evolution.c` provides SSPRK2, whole-step rollback,
actual stage/final geometry checks, positivity checks without floors,
boundary conservation ledgers, GCL diagnostics and an experimental restart
format. This format is separate from accepted production checkpoints;
restart replay qualification is still pending. The registered CFL is 0.4,
cone safety 0.9, Newton cap 8 and PCG cap 64 per Newton solve in every test.
The reported maximum CG total per RHS can exceed 64 because it sums
several Newton solves; history also sums the two RK stages.

Pressure recovery is the declared FV mean-state approximation
`P=(gamma-1)*(E-|p|^2/(2m))/V`. It is not exact subcell thermal recovery.
No extra spin state, energy reservoir, kinetic deconvolution, generator
steering, RHS mean subtraction or dynamic weight recentering is added.

## Geometry defect isolated and corrected in an opt-in path

The historical `getwfrac` clips a mathematical power-plane fraction to
[0.05,0.95]. The new moving-face formula differentiates the actual power
plane, so this clip is inconsistent with that formula once it activates.
The initial coupled implementation also encountered the legacy 2D
intersection inset and finite tiny-edge removal. Those failed sources,
accepted states and logs are retained in `fv_native_coupling_v1`.

New explicit `Voro2D_FindVCPowerStrict` and `Voro3D_FindVCPowerStrict`
entry points use the raw fraction and support a nonempty cell whose
generator lies outside it. The 2D path uses long-double plane/intersection
arithmetic without a finite inset or finite edge cutoff. Candidate-radius
certification includes weight differences; negative-offset planes postpone
the old radial early exit. Invalid/empty geometry fails explicitly.
Historical entry points retain their previous fraction policy. A native
W=0 regression compiled from parent and current sources produces identical
output bytes for the existing 1D/2D/3D geometry audit.

The independent exact fixture has x faces [0,0.1,0.5,0.9,1], raw fractions
-0.1 and 1.1, and checks actual native geometry in all three dimensions.
Maximum relative volume error is 8.88e-16, centroid error/h 4.44e-16,
degree-four monomial error/(V*h^degree) 8.33e-16, and weight-gauge error 0.
The original raw finite-difference gate **failed**, with rate errors
1.129e-6, 5.646e-7 and 2.823e-7 at fixed halving increments. The fixture
contains zero-area face births, consistent with a C1 volume and a centered
derivative error of first order in the increment. Two preregistered linear
Richardson pairs give maximum error 5.89e-10 against the unchanged 1e-7
limit. Both the original failure and the revised derivative check are
archived; this changed the diagnostic, not the PDE algorithm.

## Completed endpoints and exact Sod comparisons

The corrected v3 native binary completes all ten small coupling cases:
stationary/periodic translating contact in 1D N200, 2D 8^2 and 3D 4^3 to
t=0.25; Sod in 1D N200/N400, 2D 8^2 and 3D 8^3 to t=0.2. All six contact
noise checks pass the registered 1e-9 pressure/velocity limit. These
periodic translating contacts are separate controls, not the original
open-boundary moving-contact problem.

The same binary also completes 2D Sod 64^2/128^2 and 3D Sod 16^3 to t=0.2.
The exact Riemann solution is integrated separately across its waves on
the actual planar cell intervals, after checking transverse geometry and
velocity. Pressure is compared with physical mean pressure; continuum
subcell kinetic variance is not relabelled as thermal energy.

| Dimension and resolution | Density L1 | Mass-mean velocity L1 | Physical mean pressure L1 |
|---|---:|---:|---:|
| 1D 200 | 0.00287091 | 0.00528724 | 0.00294041 |
| 1D 400 | 0.00147138 | 0.00268281 | 0.00151114 |
| 2D 8^2 | 0.0322968 | 0.0765542 | 0.0420333 |
| 2D 64^2 | 0.00818823 | 0.0139560 | 0.00836024 |
| 2D 128^2 | 0.00421618 | 0.00719318 | 0.00422770 |
| 3D 8^3 | 0.0315495 | 0.0770914 | 0.0422172 |
| 3D 16^3 | 0.0232731 | 0.0497617 | 0.0276797 |

All these runs have positive saved states and no floor or energy repair.
Conservation including reflected-wall ledgers is within 3.81e-16 and the
normalized mass--weight invariant including its integrated residual within
2.64e-12 over the thirteen corrected runs. Sod refinement reduces all
three tabulated errors. A dimension-dependent wave timestep means these
are not bit-identical trajectories across dimensions.

## Limits and remaining work

Gresho, cold shear, the remaining original nine-test 1D suite under this
new operator, KH and the other 2D tests remain unqualified. Their earlier
results used a different operator and are not transferred here. The new
smooth phase is registered but has not run at this checkpoint. Gravity,
new long Kepler, full CUDA/MPI FV evolution and restart replay are pending.

SSPRK2 does not enforce exact finite-step GCL. The cumulative absolute 1D
Sod volume defects are about 9.98e-5 and 1.00e-4 at N200/N400; temporal
convergence and cold-flow impact still require measurement. The large
Sod runs incur 148/270/66 rejected attempts in 199/388/71 accepted steps,
so rejected work is a material performance concern.

Recorded evolution times are 6.20 s (2D64), 39.58 s (2D128), and 9.66 s
(3D16), each requesting eight OpenMP threads with up to two simultaneous
cases. Affinity overlap was not excluded. These are diagnostic timings,
not isolated or matched GIZMO benchmarks. The evolution timer excludes
initialization and final I/O; geometry timing includes initialization.
There is no established <=2x GIZMO result for this implementation.

## Evidence and reproducibility

Authoritative run root:
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`.
The tested build is `fv_native_build_v3`, parent commit
`05c51747262098e74a1856c987eb605769cd2b46`, binary SHA256
`4775b83a3533decf0791815c584f35da288b4a2386a3375bdff1a874c1514510`.
Every listed source hash matches the code published at this checkpoint.

`Research/FVNative_20261105` archives build/source manifests, frozen
sources, endpoint states, exact-reference assessments, and both failed
and corrected diagnostics. Raw source/log/data files are compressed
losslessly. Native executables stay in GPFS and are identified by hash.
The archive README specifies replay commands and outstanding scope.
