# Net-mass ALE transport and the first-RHS thermal representation

This bundled stage follows the midpoint-pressure rejection in `2e43dc7`.
Astra selected a small transport/representation test rather than a production
solver, large weight solve, or new PDE trajectory. The finite test is complete.
Conservative net-mass ALE transport closes the tested balance laws, but the
scalar mean-kinetic thermal interpretation remains unqualified. The complete
physical-accuracy, dimensional-suite, and matched-speed objective stays open.

## One specified physical transport model

For density rho, physical velocity u, face normal n, and actual geometric
normal velocity w, use the complete ALE fluxes

    Fm = rho (u_n-w_n).
    Fp = Fm u + P n.
    FE = Fm [e+|u|^2/2] + P u_n.

The energy pressure work uses physical u_n. The moving-control-volume form
follows the [AREPO finite-volume formulation](https://arxiv.org/abs/0901.4107).
No face mass flux is deleted. Fixed cell mass requires the sum of face mass
fluxes to vanish, rather than requiring every face to be material.

The test fixes uniform positive density. A small reference solve selects
weight rate lambda from

    w_n = w_n(qdot) + (lambda_i-lambda_j)/(2 distance_ij).
    L_rho lambda = b.
    (L_rho)_edge = rho A/(2 distance_ij).

The raw RHS is never mean subtracted. A pinned direct solve fixes the constant
null mode, followed by an arithmetic-mean-zero solution gauge. The solve is
an existence/reference calculation on at most 27 cells. No production weight
solver or performance result is implied. Density discontinuities and a
lambda-dependent Riemann state remain outside the test.

The cell-mean velocity and generator motion are distinct physical quantities.
For uniform-density affine velocity with exact gradient G and cell centroid c,

    ubar_i = p_i/m_i.
    qdot_i = ubar_i + G(q_i-c_i) = u(q_i).

This equality is checked using the known analytic affine field. No general
reconstruction from arbitrary states is selected. The distinction introduces
no mesh-regularization velocity or prescribed orbit.

## Native geometry and integration accuracy

The original native interval/Voro/strict-3D routines replay the previously used
1D 16-cell, 2D 9-cell, and 3D 27-cell geometries. Positions and volumes match the
existing geometries. A compile-time-only diagnostic hook records actual face
vertices and positive surface quadrature. There are 32, 108, and 12,348 surface
points, respectively. No source from the native tessellation library changes.

Two-point line Gauss quadrature and a positive 3x3 Duffy triangle map integrate
degree-three uniform-density affine energy fluxes. The triangle map is the
existing `gfs_dg_quad` construction generalized to a 3D face. No DG state,
evolved moment, or density closure is enabled. Higher-order face integration
is also discussed for cold shear in
[Zier & Springel](https://arxiv.org/html/2205.07916#S4.SS2). The citation does not
qualify the present model or establish a runtime for LagEunha.

Quadrature reproduces native face area and centroid. Divergence-theorem
integrals give cell volume, centroid, and central second geometric moment Q.
The latter is used only as a diagnostic integral. With the optional hook
disabled, original and modified geometry-audit source compile to bit-identical
binaries and produce bit-identical 1D/2D/3D diagnostic output.

## Ten frozen controls close the transport balances

The controls are translation and rest linear pressure in every dimension,
plus rigid rotation and stationary affine shear in 2D/3D. No nonzero shear or
rotation control is claimed in 1D. Gamma matches the stored native value.
Uniform-pressure cases use rho=1 and P=1e-8. The linear-pressure case uses a
positive field with base pressure two and gradient 0.2 along x.

The physical affine field is continued onto the owned-cell patch. The geometry
keeps its original fixed periodic image lattice. Nonperiodic field values at
the cut enter explicit external mass/momentum/energy boundary ledgers. The
control does not implement a production shear-periodic boundary condition or
claim conservation without the external cut terms.

The reference solves retain nonzero surface mass transport during rotation
and shear while each cell has zero net mass rate to roundoff. The centroid
transport identity is independently checked as

    m_i cdot_i = p_i - sum_faces integral x Fm.

Euler momentum rates use the exact initial temporal derivative of the analytic
field. A rigid-rotation velocity field with constant pressure is not claimed
to be a force-free stationary Euler solution. For incompressible affine data,
the initial Euler derivatives include Gdot=-G^2 and the corresponding uniform
velocity derivative. Galilean boosts transform all flux components and the
geometric face speed together.

The largest algebra/balance defect among all ten controls is 5.6e-16.
Galilean flux/thermal-rate defects are at most 7.4e-16. Native position/weight
geometry replays at the two fixed increments 1e-5 and 5e-6 independently check
volume, centroid, and central-moment derivatives. The largest volume/centroid
derivative defect is 1.24e-9. The largest relative kinetic-variance derivative
defect is 4.72e-8. No extra increment search or RHS compatibility adjustment
was used.

## Initial energy choice does not remove the first-RHS issue

The old Hopkins exact-energy initialization experiment is not repeated. The
new test asks whether choosing initial discrete energy `U_physical+Kmean`
removes the problem from the first update. No actual checkpoint is converted.
Evaluate the physical ALE energy flux and define its mean-kinetic proxy rate

    Uproxy_dot = Edot - ubar dot pdot + |ubar|^2 mdot/2.

The final term vanishes for exact fixed mass and is retained numerically.
For an actual cell-integral interpretation, the unresolved kinetic variance is

    Kvar = (rho/2) trace(G Q G^T).

The diagnostic physical thermal rate is zero for the registered incompressible
controls. The measured proxy rate satisfies

    Uproxy_dot = Kvar_dot.

The maximum discrepancy is 1.55e-17. Kvar_dot is checked independently from the
native perturbed geometric moments and the physical affine gradient derivative.
It is not assigned to a new reservoir or subtracted from an evolved energy.

| Dimension | Field | Absolute Uproxy-rate L1 / volume | Signed proxy-rate sum | Absolute face mass transport |
|---|---|---:|---:|---:|
| 2 | rotation | 6.11481e-4 | +3.14901e-4 | 1.18838 |
| 2 | stationary shear | 1.46031e-3 | -3.51186e-4 | 1.27011 |
| 3 | rotation | 8.07023e-4 | +1.29710e-4 | 1.07286 |
| 3 | stationary shear | 1.24334e-3 | -7.30202e-5 | 1.12607 |

Translation and rest linear pressure have zero proxy rate to roundoff.
The nonzero rates have both signs and must not be described as positive
dissipative heating or measured entropy production. The surface transport
column is an absolute sum; the cell net mass rates remain at roundoff.

One half-length scaling of the same local shear pattern reduces the proxy
rate per volume by a factor 0.25 in both 2D and 3D. The result identifies an
h^2 local spatial-representation term. It is not a full fixed-domain PDE
convergence result and does not meet a cold-flow error budget automatically.
No value is extrapolated to the N128 Gresho or previous Kepler trajectories.

## Current decision and open closure

The tested ALE transport identities pass. A model interpreting
`E-|p|^2/(2m)` directly as physical thermal energy is not adopted. Fixing only
the initial projection leaves the measured variance-rate term in the first
physical flux. The result does not prove that all finite-volume cold-flow
methods fail or that every subcell kinetic contribution is irreversible heat.

A next model must specify the conserved state, velocity location,
reconstruction, kinetic-energy interpretation, and pressure recovery together.
A literal cell-average model and a particle point-quadrature model cannot be
interchanged by renaming pressure after a failed test. No additional evolved
moment, thermal reset, limiter, pressure-implicit method, N128 cost/spectral
test, or new production trajectory is introduced by the present bundle.

Authoritative evidence is under
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/ale_balance_feasibility_v1`.
Compact source, native quadrature, reference data, finite differences, and
checks are archived in `Research/ALERepresentation_20261105`. The original
fluid executable and manuscript benchmark remain unchanged.
