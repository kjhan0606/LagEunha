# Voronoi geometric and material moment rates in Q2

This diagnostic addresses the unresolved W=0 fixed-mass affine/spin closure.
For each saved Q2 inner-disk polygon, assume uniform density and define the
central mass second-moment tensor `S` and mass centroid `C`. There are two
different instantaneous derivatives:

* If the polygon remains uniform density while its generators move and its
  mass stays fixed, its **geometric** derivative follows its moving faces:
  `Sdot_geom=(m/V) sum integral[((x-C)(x-C)^T-S/m) w_n] ds`.
* For a material region with no mass crossing and affine velocity gradient
  `G`, its **material** derivative is `Sdot_material=G S + S G^T`.

The first expression is evaluated with the actual stored generator velocities,
the exact Voronoi face speed and two Gauss points per face. Independent
centered differences of polygon moments after moving **all** generators
verify the face result: relative Frobenius RMS error is 4.92e-8 at t30 and
5.88e-9 at t47.5 for the smaller tested step. The stored-gradient result is:

| Time | `mean(||Sdot_geom||/(Omega tr S))` | `mean(||Sdot_material||/(Omega tr S))` | `mean(||difference||/(Omega tr S))` |
| --- | ---: | ---: | ---: |
| 30 | 0.8710 | 0.8424 | **1.7132** |
| 47.5 | 0.8699 | 0.8426 | **1.7118** |

Means are mass-weighted over the same 3,931 cells initially in `.5<=r<.75`;
matrix norm is Frobenius and `Omega` uses the instantaneous softened circular
frequency. These are instantaneous rates, not accumulated error or heating.
The difference is comparable to the sum of the two rate magnitudes. Replacing
the saved gradient with the analytic softened Kepler gradient changes the
mean difference only to 1.7133 and 1.7119, respectively. The saved versus
analytic gradient RMS discrepancy is 5.51e-4 and 1.50e-3 relative at the
two checkpoints. The latter control still uses stored generator velocities
for geometric face motion, so it is not a fully analytic evolving solution.

Under these explicit assumptions, **identifying a transported material
second moment with the instantaneous uniform-density Voronoi polygon moment
does not give the same evolution** in this shear flow. The difference is a
constraint on the proposed model: it must represent the material/geometry
relationship or derive a consistent non-Euler discrete closure. It does not
prove that the existing native fixed-mass variational solver is invalid, that
spin cannot work, or that all nonuniform-density closures fail. Actual native
gradients can correspond to a nearby RK stage, the comparison is instantaneous,
and the inner disk is bounded by other cells rather than an independent
domain. There is no viscosity, steering or new simulation in this analysis.

See [stored-gradient receipt](KEPLER_Q2_MATERIAL_MOMENT_GAP_2026-10-29.json),
[analytic-gradient receipt](KEPLER_Q2_MATERIAL_MOMENT_GAP_ANALYTIC_2026-10-29.json),
[stored-gradient reader](audit_kepler_q2_material_moment_gap_20261029.py),
and [analytic-gradient reader](audit_kepler_q2_material_moment_gap_analytic_20261029.py).

## Can instantaneous weight rates match all three fixed-mass moments?

At each saved W=0 state, form the linear response of the uniform-density,
**fixed-cell-mass geometric** `Sdot` to signed power-weight rates. A face
between sites i and j has weight-induced normal speed
`(Wdot_i-Wdot_j)/(2 distance)`. Integrating the three central-moment
components over each face yields a sparse Jacobian; all 3,931 inner cells
contribute three equations, while rates for all 4,646 sites touching those
cells—including neighbors outside the selected region—are unconstrained
variables. The least-squares norm uses mass-weighted Frobenius components,
normalized by `Omega tr(S)`. No Wdot sign, magnitude, pressure-work or
weight-feasibility restrictions are imposed. A random Jacobian-vector probe
against independent power-cell finite differences has relative RMS error
1.46e-6 at t30 and 1.24e-8 at t47.5.

The initial 1,500-iteration runs hit the iteration limit. Increasing the
limit produced LSMR convergence status 2 after 5,141 and 5,795 iterations.
The relative **weighted RMS residual** `||J Wdot - (Sdot_material -
Sdot_geom)|| / ||Sdot_material - Sdot_geom||` is **0.99799 at t30** and
**0.99600 at t47.5**. The relative normal-equation residuals are 2.63e-7
and 1.17e-7. Thus this specific unconstrained numerical fit removed less than
0.5% of the mismatch norm even though its weight rates were free on the
outer neighbors. It is an instantaneous linearized numerical result, not a
rigorous singular-value proof, a complete nonlinear feasibility analysis,
or a long-time evolution.

This calculation deliberately held each cell's mass fixed while changing its
weight. The user also permits weight-linked **mass exchange** for nonzero
changing weights; such a model transports mass and its moments and requires
different coupled equations. Therefore this result does **not** rule out
evolving weights, nonuniform density, or the approved variable-mass path; it
only argues against using scalar weight rates alone to make all three
fixed-mass uniform-polygon moments follow the chosen affine material rates in
these Q2 states. No native simulation, generator steering or periodic boundary
was introduced. See [converged receipt](KEPLER_Q2_WEIGHT_MOMENT_FIT_LONG_2026-10-29.json),
[short receipt](KEPLER_Q2_WEIGHT_MOMENT_FIT_2026-10-29.json), and
[long reader](audit_kepler_q2_weight_moment_fit_long_20261029.py).
