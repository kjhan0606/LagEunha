# Material volume transport and a constrained-weight candidate

This continues `70a2694` in the same LagEunha repository. The full objective is
still open. The earlier 1D contact initialization defect is resolved separately;
the present investigation concerns multidimensional volume work and vortex
accuracy. Production pressure, viscosity, weights, and time integration have
not been replaced. New native hooks are read-only frozen-state diagnostics.

## Actual geometry identifies a transport error

For the current unweighted operator, define the local affine transport moment

    K_i = sum_faces (bj tensor r) / V_i.

Its volume rate for an affine velocity gradient G is `Vdot/V = K:G`, while
physical divergence is `trace(G)`. Actual geometry gives a linear-pressure
moment `sum(bi tensor r)/V = I`, symmetric K, and `trace(K)=D`. These imply
correct linear pressure, rigid rotation, and dilation patches. They do not
imply `K=I` or correct local volume transport in general incompressible strain.

A repeated irregular site pattern gives:

| Dimension | Cells, coarse/fine | Incompressible affine volume-rate L1, coarse/fine |
|---|---|---|
| 1 | 16 / 32 | roundoff / roundoff |
| 2 | 9 / 36 | 0.0810341 / 0.0810341 |
| 3 | 27 / 216 | 0.0663607 / 0.0663607 |

The 2D/3D pattern is exactly repeated under refinement. Signed total volume
rates cancel to roundoff, as do direct affine viscous heat and fitted physical
divergence. The effect therefore precedes direct AV heating and is not a failure
of global energy/volume bookkeeping. This is a local affine consistency test,
not a theorem that all Gresho solutions fail to converge.

Ten frozen native Gresho states, including the timestep controls, separate the
same effect in the actual solver. At 128² t=0.5, exact Gresho generator velocities
give volume-rate L1 0.54168 versus fitted divergence 0.03575. Contracting the
exact analytic gradient with K already gives 0.52641. In the smooth r=0.25--0.35
annulus, the corresponding values are 1.43499, 0.03958, and 1.43140. Exact affine
core rotation remains at roundoff. The effect persists in the dt/4 control.

Force error is split linearly into the error with exact **nodal** pressure and
the force caused by the evolved pressure deviation. These are force-density
norms; the pressure oracle is not asserted to be a cell-average projection.

| State | Exact-nodal-pressure discretization error | Evolved pressure-noise contribution |
|---|---:|---:|
| 64² t=0.5 | 0.03560 | 3.79565 |
| 128² t=0.5 | 0.01569 | 4.92320 |
| 64² dt/4, t=1 | 0.04109 | 1.41659 |

CPU/GPU force differences are below 3e-17, thermal RHS agreement with `-P Bv`
is below 6e-17, and linear-pressure moment errors are below 1.3e-15. The main
remaining force disturbance is the evolved pressure field rather than an
incorrect implementation of the variational pressure gradient.

## Exact-flow geometry control

Birth generators are also mapped by the exact stationary Gresho circular flow,
with no force integration or PDE step. The original native tessellator measures
the new volumes. Keeping each birth mass, the geometric density L1 from one is:

| Mapping time | 64² | 128² | 256² |
|---:|---:|---:|---:|
| 0.25 | 0.01248 | 0.00493 | 0.00204 |
| 0.5 | 0.07205 | 0.04216 | 0.01918 |
| 1 | 0.13442 | 0.11888 | 0.09637 |

This control isolates geometric mismatch without numerical orbital errors or
AV. The native input retains its birth thermal data solely for geometry
measurement; its force and heat outputs are not treated as a physical state.
The separately calculated fixed-entropy prediction is
`P/P_birth=(V_birth/V_actual)^gamma`. At mapping t=1 its relative L1 is
0.2392/0.2114/0.1712. These are predictions along a prescribed geometry path,
not evolved hydrodynamic results or an energy-conservation test.

The old material-face proof already rules out exact pointwise material faces
from scalar weight rates under general shear. It does not rule out net-volume
correction. That old proof and the failed higher-moment branches were consulted
before selecting the next finite test; they were not rerun.

## Candidate selected by the bundled Astra review

Let `L=dV/dW`, with interior coefficients `A/(2d)` and a constant null mode.
For orthonormal traceless symmetric basis tensors S, solve geometric fields:

    L psi^b = V (I-K):S^b.

There are two fields in 2D and five in 3D. The velocity action is local:

    (Rv)_i = sum_b psi_i^b (G_i:S^b) + (2/D) W_hat_i trace(G_i).

Both psi and W_hat use arithmetic-mean-zero gauges. The candidate must couple
all of the following, rather than change the thermal divergence alone:

    Wdot = Rv.
    C = B + L R.
    m vdot = C^T P.
    Udot = -P Cv.

The coefficient construction is global, but R and `R^T(LP)` act on local fit
stencils. A common shift of psi is physically significant for varying strain,
so fixing its gauge is part of the definition. At W=0 the active constraint can
still have R nonzero; full return to the old operator requires disabling R.
Mass remains fixed in this proposed constraint model. No canonical Hamiltonian,
symplectic KDK, nonlinear Euler accuracy, or new physical trajectory is proved.

## Completed necessary gates

1. Six small geometries in 1D/2D/3D pass compatibility, gauge, constant pressure,
   translation, rotation/dilation/general affine transport, force/work adjoint,
   and momentum checks. Affine transport errors are at most 7.1e-15. Interior
   linear-pressure tests have four 2D and eight 3D samples; periodic-cut work is
   recorded explicitly. Periodic owner angular momentum is not called an invariant.
   In 1D at W=0, R is identically zero.
   Independent frame changes of the exported 2D/3D geometry preserve R, C and
   pressure force to 1.7e-16; Galilean changes in R are below 3e-18.
2. For a fixed physical wavelength periodic shear on 6²/12²/24², the original
   volume-rate/fit divergence gap stays near 0.32. The candidate gives
   0.08264/0.03890/0.01877. Physical divergence L1 is
   0.12773/0.06458/0.03251. The varying-strain commutator is included; its exact
   decomposition residual is below 7e-17. This is one necessary spatial test.
3. Coupled position/weight finite differences in 1D intervals and actual native
   2D/strict-3D geometry pass at W=0 and at 0.02 h² weights. Best relative errors
   range from 1.4e-13 to 1.7e-10. A common weight shift changes volume by at most
   1.1e-14 and the constructed operator by at most 6.5e-15. No clipping or jitter
   was added. A first verifier incorrectly required further improvement after
   an exactly linear 1D result already reached roundoff; the unchanged accuracy
   threshold is retained, with that verification defect explicitly recorded.

Dense NumPy eigensolves provide independent small-matrix reference evidence.
They do not evolve a fluid and are not the production coefficient solver or a
performance measurement. The geometry comes from the original C routines.

## First native cost screen: not qualified

The saved 128² t=0.5 Gresho geometry has 16,384 cells and 98,306 directed faces.
A read-only original-native export supplies exact owned/image aliases and
relative positions. Every face is retained. Two earlier attempts to infer
padding aliases from almost-zero edge directions failed input conversion and
never executed the cost solver; those technical records are preserved.

The bounded C/OpenMP Jacobi-PCG screen uses actual eight threads, tolerance
1e-10, maximum 256 iterations per field, and a total extra-work budget of one
median original full native step: 0.076127 seconds. The first field stops after
175 iterations at 0.076634 seconds with residual 5.63e-5. The second field,
R/R-transpose evaluation, weighted geometry, and a new PDE trajectory are not
released by this result. Iterations were not increased after the failure.

This implementation fails the finite cost/convergence screen. It does not
disprove the geometric model. Its inherited `OMP_WAIT_POLICY=passive` and many
short parallel regions motivated one further bounded execution test selected
by Astra. The OpenMP setting is an implementation hint, so its performance
effect was measured rather than assumed. See the
[OpenMP wait-policy definition](https://www.openmp.org/spec-html/5.1/openmpse64.html).

## OpenMP crossover and fair native baseline

One instrumented C binary was tested in three PASSIVE/ACTIVE crossovers with
the same graph, eight cores, binding, zero initial guess, numerical tolerance,
256-iteration ceiling, and 0.076127-second total budget. True residuals were
recomputed after every attempt, including truncated solves.

| Policy | First-field iterations | First-field residual | Median iteration time |
|---|---|---:|---:|
| PASSIVE | 180, 180, 186 | 2.7e-5--3.9e-5 | 0.411 ms |
| ACTIVE | 256, 256, 256 | 2.6414e-7 in all runs | 0.112 ms |

ACTIVE reduces execution overhead, but all runs hit the unchanged iteration
ceiling without meeting 1e-10. Under the preregistered stage rule, no persistent
team implementation, iteration increase, or new preconditioner search followed.
The second field, R/R-transpose, weighted geometry cost, and candidate PDE remain
unqualified. The current CPU coefficient implementation is held.

The original immutable native executable also continued the same accepted
t=0.5 checkpoint to exactly t=0.51 under both policies, with the same eight-core
settings and one final output. Both runs took 14 accepted steps. All final
x/y/m/v/W/weight-momentum/U/V fields are bit-identical; per-cell masses remain
fixed. Results and relative energy change from the common start are identical.

| Policy | Native phase sum | Median step after the first | Total process wall time |
|---|---:|---:|---:|
| PASSIVE | 1.16294 s | 0.065530 s | 2.07201 s |
| ACTIVE | 0.78282 s | 0.043839 s | 1.67569 s |

This short matched control measures 1.486x native-phase and 1.237x process-wall
speedup. It is not a long-run or GIZMO throughput result. The first incomplete
coefficient attempt alone costs 0.694 of an ACTIVE native median step, before
the second field or later work. The faster baseline therefore does not justify
comparing the candidate only to the older, slower timing.

Future controlled native tests should explicitly use ACTIVE with reserved CPU
resources; historical run environments remain unchanged. This measured runtime
improvement is separate from the unresolved physical volume-transport problem.

## Independent frozen-graph reference, analysis only

An independent sparse factorization of the same 128² t=0.5 native graph supplies
mean-zero reference coefficients for the planned residual-to-physical-error
analysis. The two relative reference residuals are 4.6251e-14 and 5.8694e-14.
The maximum entry of `L-L^T` is 4.44e-15. Saved matrix, coefficient, and RHS
arrays reproduce the reported residuals in an independent post-run check.

Ritz values for the symmetric part of the matrix suggest positive-spectrum
condition ratios of about 12,596 without scaling and 4,097 with Jacobi scaling.
The measured Ritz residuals accompany the estimates. The estimates are not
certified global spectral bounds and cannot certify a physical error budget.
The factorization timing is reference-only and is not a production cost result.
The script performs no fluid evolution and accepts no PDE steps.

Compact results and SHA256 provenance are in
`Research/MaterialVolume_20261105/material_volume_reference_v1`. The preceding
256-iteration CPU cost failure remains unchanged. The planned finite stage
propagates coefficient errors into volume work, force, and face displacement
without increasing the iteration cap or relaxing the tolerance. The completed
follow-up is now in `UNIFIED_HYDRO_COEFFICIENT_ERROR_2026-11-05.md`. The follow-up
also identifies a frozen velocity/pressure-noise amplification with accurate
reference coefficients. No candidate PDE or manuscript result is released.

## Provenance and remaining work

Run/build root: `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`.
Compact evidence is in `Research/MaterialVolume_20261105`; large native states,
raw CSVs and geometry files remain on GPFS with hashes. The native diagnostic
build is a partial relink. No old accepted state or ongoing job was replaced.
The manuscript retains its previous verified baseline.

The next stage requires a method/implementation decision after this completed
bounded cost failure. Do not resume the same CPU solve with a larger cap or
silently relax its tolerance. Full weighted geometry cost, constrained time integration,
physical vortex/shock accuracy, dimensional PDE regressions, remaining 2D tests,
new-method Kepler, and matched GIZMO throughput remain open.
