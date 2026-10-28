# Kepler disk: fixed-geometry Laguerre volume feasibility, 2026-10-26

## Why this probe was run

The C010 fixed-mass, zero-weight Voronoi run shows about 10.49% mass-weighted
absolute density change in its first disk ring (`2 <= r < 2.25`) at `t=6`.
An independent reconstruction in which every generator follows its **exact**
softened circular orbit reproduces the signed per-cell density change with
correlation `0.99999987`. The partition is not material under Kepler shear.
Changing the mesh weights is therefore a concrete way to test whether the
geometric part of this error can be removed without steering generators.

## Mathematical test

For each generator (i), define the power cell by minimizing
`|x-x_i|^2-w_i`. Keep `w_i=0` outside the selected active region. At the exact
ballistic position of every generator, solve for active weights so that every
active cell's area equals its `t=0` native Voronoi area. Along a finite shared
face (ij), the derivative of cell area with respect to weights is
`L_ij/(2|x_j-x_i|) * (delta_w_i-delta_w_j)`. The assembled matrix is a sparse
weighted graph Laplacian with zero-weight exterior cells as Dirichlet boundary.

The independent reference implementation builds the regular triangulation as
the lower convex hull of lifted points `(x,y,x^2+y^2-w)`, computes dual power
vertices, polygon areas, and each current Jacobian, then performs damped Newton
steps with Jacobi-preconditioned conjugate gradients. It validates its
unweighted power-cell areas against SciPy's separate planar Voronoi routine:
median relative difference `4.04e-13`, 99th percentile `3.23e-12`.

## Results

All figures below use the **same 71,080 generators** and the same initial
volume targets. The disk-only solve has 21,436 active cells with
`2 <= initial r < 8`; the inner taper and exterior remain at zero weight.

| Exact-orbit time | Starting weights | First disk-ring median absolute relative volume error | Final maximum relative error among all active disk cells | Outcome |
|---|---|---:|---:|---|
| 6 | zero | 0.0798868 | 1.48e-9 | one linear solve plus two full Newton steps |
| 6.250024729493809 | converged t=6 weights | 0.0175685 | 4.85e-12 | three full Newton steps |

No active disk cell disappears or has a nonpositive area in either final
diagram. This positive-area result does **not** by itself meet the native
weight-difference guard described below. The zero-weight reference error in
the first ring at `t=6.25` is
`0.107883` by this median metric. The absolute weight change between the two
solved times has median `1.02e-5`, 95th percentile `1.04e-4`, and maximum
`5.92e-4`; its L2 norm is 10.2% of the t=6 weight L2 norm. This 0.25-time
step spacing is far larger than the live RK4 step (~2.5e-5).

The one-time CPU reference calculations took 12.9 s for t=6 and 15.0 s for
t=6.25, including 71k-site power-diagram construction. CG used roughly
300--420 iterations per Newton step. Those are **analysis costs**, not an
acceptable hydro evolution cost or a GPU benchmark. A per-stage implementation
must warm-start and demonstrate near-AREPO throughput before adoption.

### Inner-taper caveat

Extending the active target to `1 <= r < 8` adds 1,076 inner-taper cells. The
linear solve reduces disk-ring median error from 0.07989 to 0.00147, but the
current nonlinear iteration stalls with a maximum active-cell relative error
of 0.254, located at initial `r~1.09`, `1.38`, and `1.58`. A finite-difference
Jacobian check has 0.709 relative L2 discrepancy over this expanded region,
whereas the disk-ring part agrees at `1.6e-6`. This does **not** prove that a
solution for the taper does not exist. Near-degenerate geometry, active-set
changes, the zero-weight exterior condition, and the reference implementation
need separate investigation. The disk-only result is much cleaner and is the
appropriate first native implementation target.

### Native face-validity guard

The existing `GFS_DYN_W2` code uses `|w_i-w_j|/d_ij^2 < 0.8` as its default
face guard. The exact-volume disk solution exceeds that guard on only four
symmetry-related faces out of about 65,000 faces touching active cells, but
all four lie across the sharp `r=2` boundary (`r_i~1.990`, `r_j~2.016`). At
`t=6`, their maximum ratio is `0.9945`; at `t=6.25`, it is `1.2179`.
The weighted polygons still have positive area, but the native guard prevents
using this exact solution unchanged.

Extending the active region to `r>=1.9` also converges geometrically at
`t=6.25` (maximum relative volume error `3.24e-8`), yet it spreads the guard
issue: 20 faces exceed 0.8 and four exceed 1. A simple linear radial taper
of the `r>=2` weights over width 0.025 eliminates the guard violations, but
the first-ring density-error 95th percentile rises to 6.4% and the largest
volume error reaches 20.9%. Directly reducing only the four offending cell
weights also meets the guard but creates up to 26.2% density error in those
cells. Wider ramps perform worse in the first disk ring. Neither ad hoc
alteration is suitable for promotion as a physical solver.

The evidence points to a constrained volume/face-motion problem at the disk
boundary, including conservative mass flux where a face cannot follow the
local fluid. It must be formulated and tested in the native flux operator; a
static weight cap would hide the local error rather than remove it.

### Native checkpoint cross-check and volume-rate probe

Repeating the disk-only solve on the **actual** C010 generator positions from
the audited t=6 and t=6.250024729493809 checkpoints gives the same result:
maximum active-cell relative target-volume errors `1.48e-9` and `4.85e-12`.
It identifies exactly the same four boundary faces over the guard, with
maximum `|dW|/d^2` ratios `0.9950` and `1.2188`. Thus the boundary issue is
present in the native run and is not an artifact of exact-orbit reconstruction.

At the actual t=6.25 state, moving generators briefly with their saved
velocities while holding weights fixed gives a reference geometric volume
rate `b_i`. Solving the current power-face graph Laplacian
`J dot(w) = -b` makes the disk volume stationary to first order. In a separate
forward Euler geometry check with `delta t=0.001`, the first disk-ring median
absolute relative area change is `8.27e-5` if weights stay fixed, `2.99e-7`
with 32 CG iterations, and `7.23e-8` with the fully converged rate. The
95th-percentile changes are `2.72e-4`, `1.32e-6`, and `1.73e-7`. Sixteen CG
iterations give a median `1.45e-6` and 95th percentile `6.08e-6`.

Cold convergence to relative linear residual `9.5e-11` takes 322 CG
iterations; a secant rate from t=6 and t=6.25 still takes 308. The 32-step
linear residual is `0.00918`. These results suggest a truncated GPU solve
could be useful, but there is no native RK4, conservation, or GPU timing
result yet. The four guard violations remain even if the volume rate is
accurate.

Differentiating the power bisector gives the **normal geometric face speed**
at face centroid `c_f`, with `d=|x_j-x_i|`:

`v_face,n = [x_j·v_j - x_i·v_i - (v_j-v_i)·c_f + (dot(w_i)-dot(w_j))/2]/d`.

Summing each exact face length times this speed reproduces the independent
centered finite-difference cell-volume rate at relative L2 error `2.50e-7`.
The solved rate makes the summed face-volume rate negligible at the same
finite-difference accuracy. The existing GPU `GFS_W2_LINEAR` rate term plus
`GFS_LAGUERRE_ROTATION` centroid correction algebraically implements this
normal speed, provided the geometric face velocity is actually selected for
the work term. This validates the geometric conservation-law formula, not the
full fluid conservation update.

### Signed weights and native representation

At the actual t=6.25 state, 83.5% of active disk weights in the zero-exterior
gauge are negative (`min=-0.0036761`, `max=0.0049976`). A common offset of
`+0.0036761` to **all** 71,080 weights makes them nonnegative without changing
the power diagram: the independently computed cell areas agree within
`4.85e-12` relative. This is a gauge transformation, not a change in physics.
The shifted range is `0` to `0.0086737`.

The native GPU path still cannot accept these values unchanged. Its current
`kappa=0` gate clears weights, `GFS_W2_LINEAR` only activates for positive
hydro-derived kappa and limits weights, and a current face-containment ceiling
is exceeded by exactly the four disk-boundary sites at `r~2.016`. The common
offset does not cure that gauge-invariant containment failure. A new native
mode must preserve a consistent common gauge, skip unrelated hydro-weight
clipping, and resolve the four boundary faces physically before a Kepler run.

As an early stencil-cost check, weighted power neighbors touching the active
disk differ little from zero-weight neighbors at t=6.25: 80 versus 72 faces
span more than two linked cells in one direction, and none span more than
four. This supports a bounded native tessellation test but is not a GPU
performance or candidate-completeness proof.

## What has and has not been established

For these two fixed exact-orbit point configurations, positive-area Laguerre
cells with essentially the initial volume exist for every active **disk**
generator under the stated zero-weight boundary. This is strong geometric
evidence for time-dependent disk weights. It does **not** establish a
continuous-time solution over t=6--10, stable hydro integration, correct
weighted face velocity, mass/momentum/energy conservation, or speed.

The next native-code gate should derive `dot(w)` from the cell-volume rate at
each RK stage and use the same weighted face motion in the flux and energy
operators. Keep generator trajectories unsteered. Preserve the reduction rule
`w=0 => Voronoi with fixed cell mass`. If nonzero weights permit face mass
exchange, track it exactly in conservative variables and tie any relaxation of
fixed mass to the weight evolution. Benchmark the added linear solve and mesh
rebuild cost on the GPU before a long Kepler run.

Existing native paths need careful reuse: `SEDOV_LAGVOL` does one limited
volume corrector and replaces density with a target-volume value, so it does
not certify agreement with the realized mesh. `GFS_DYN_W2` solves a different
per-face HLL speed fitting problem and is disabled on GPU. `GFS_W2_LINEAR`
already transports a weight rate into the GPU face-velocity/work operator but
currently requires positive hydro-derived kappa; its coupling is a useful
starting point for a conservative GPU implementation.

## Reproduction

- `laguerre_linear_volume_feasibility_t6_20261026.py` and JSON: frozen-topology
  first-order solve for the wider `1 <= r < 8` region.
- `laguerre_nonlinear_volume_probe_t6_20261026.py` and JSON: independent
  power-cell geometry at zero and first-order weights.
- `laguerre_newton_volume_t6_20261026.py`: nonlinear geometry solve.
- `LAGUERRE_NEWTON_VOLUME_R2_T6_2026-10-26.{json,npz}` and
  `LAGUERRE_NEWTON_VOLUME_R2_T625_2026-10-26.{json,npz}`: reproducible results
  and warm-start weights. Exact source snapshot is
  `/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core_refinement/fine_integrated/khout.000001.dat`.
- `laguerre_weight_validity_t6_t625_20261026.py` and
  `LAGUERRE_WEIGHT_VALIDITY_FACES_2026-10-26.json`: power-face guard ratios and
  the four boundary faces that exceed the native guard.
- `laguerre_boundary_taper_probe_t625_20261026.py` and
  `LAGUERRE_BOUNDARY_TAPER_T625_2026-10-26.json`: static ramp screening; not
  a proposed evolution rule.
- `LAGUERRE_NEWTON_VOLUME_R2_{T6,T625}_ACTUAL_2026-10-26.{json,npz}` and
  `LAGUERRE_WEIGHT_VALIDITY_FACES_ACTUAL_2026-10-26.json`: the same geometric
  solve and guard audit using native C010 checkpoint coordinates.
- `laguerre_weight_rate_actual_t625_20261026.py` and
  `LAGUERRE_WEIGHT_RATE_ACTUAL_T625_V3_2026-10-26.json`: differential
  weight-rate, truncated-CG geometry, and face-volume-rate checks.
- `laguerre_gauge_contract_actual_t625_20261026.py` and
  `LAGUERRE_GAUGE_CONTRACT_ACTUAL_T625_V2_2026-10-26.json`: common-gauge
  equivalence, containment failures, and linked-cell neighbor spans.

The independent solver changes no live C010 setting or native source file.
