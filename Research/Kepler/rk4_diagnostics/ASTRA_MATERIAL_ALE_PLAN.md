# Astra structural decision: compatible material map plus conservative Laguerre remap

Current implementation/evidence: `MATERIAL_ALE_RESULTS.md`. Local smooth
material RHS, finite affine mass-weight common-refinement controls and the
native CUDA helper are implemented. The first bundle remains incomplete:
global shared kinematics, curved Kepler remap/compression and composed
density-wave/sliding checks are still required. No runtime/long pass.

2026-10-21. One bundled review of prior split/central/reconstructed-HLL
failures and the affine-face kinematic obstruction. Driver implements native
C/CUDA helpers and bounded native checks. The old nodal pair-flux family is
closed; no damping-factor, limiter or nullspace tuning of those candidates.

## Compatible Lagrange state

    x_h(X,t)=sum N_a(X) x_a(t), v_h=xdot_h,
    F=dx/dX, J=det(F)>0, dm=rho0(X)dX,
    M_ab=integral N_a N_b dm, K=.5 v^T M v, U=integral e dm.

e is physical specific internal energy. Use the same kinematic space and
mass matrix for all represented rotation/strain energy. Never define thermal
energy by E-Kaffine and import the unresolved initial kinetic remainder.

For symmetric Cauchy stress sigma=-P I+tau_num,

    (M vdot)_a=-integral sigma grad_x(N_a)dV
                + integral N_a g dm + boundary traction,
    edot=(sigma:grad_x(v_h))/rho.

Use identical quadrature, stress and velocity gradient for force/work/heat.
Then Kdot+Udot is boundary work plus gravity power. Coordinate/constant
reproduction, symmetric M and symmetric stress give the actual angular
identity L=sum x_a cross (Mv)_a, without an independent spin correction.

## Initialization, shock and geometry requirements

Keep the actual input velocity function in a finite reference-coordinate
kinematic enrichment along with constant/affine modes; coefficients evolve,
never reset to the benchmark. Do not delete kinetic energy or send projection
loss to heat. Match initial K/U/E with independent integration first.

Smooth tau_num=0 is not a shock solver. A full Euler Riemann-derived
compression stress on material subfaces must provide symmetric stress,
nonnegative dissipation and identical force/work heating. If this is not
implemented, report the smooth-only limitation; no Sod/Sedov/runtime pass.

Material cells are auxiliary true flow maps, not Euclidean power cells.
Accepted geometry is Laguerre. Scalar Wdot cannot remove the tangential
strain mismatch of a planar power face (`AFFINE_FACE_KINEMATICS_RESULTS.md`).
Transport between material and target geometry is an Euler remap, not a
hidden mismatch of thermodynamic and geometric volume.

## First bundle includes remap

Use the original24 actual polygons: rigid/affine controls and exact original
Kepler+same uniform velocity perturbation. Measure thermal chain, actual map
volume derivative, stress torque, gravity work and boundary traction/work.
Include one finite general shear deformation, material-to-Laguerre common
refinement, and a return-remap control. Transfer m,p,E and actual angular
integrals, retaining linear momentum spatial information. Distinguish overlap
integration, kinetic projection, thermal and conservation errors.

Common-refinement field fragments may initially remain uncompressed so these
errors can be measured separately. That is not a scalable long-time state:
compression/remap must later control cold error. Never make E appear conserved
by automatically transferring lost kinetic projection energy to U.

Target weights must couple the actual remapped mass to

    m_i+mu_i W_i=m_i0.

The screened rate solve is reusable reasoning but is not a finite-remap
existence certificate. Do not hold W while dropping actual transfer.

After this, test density-wave response with a fixed-target remap, not only a
stationary material-frame state, and screen the historical cs/h sliding mode.
No full71080-cell evolution until these contracts, shock handling, actual
geometry, weight coupling and native MPI/runtime/restart all pass.

Stop on material/power-cell identification, force/work-volume mismatch,
unexplained angular source, hidden remap kinetic-to-heat exchange, sliver or
geometry loss, sliding growth, or deletion of initial kinetic energy.
This distinguishes the proposal from old Voronoi VPH force and from the old
Delaunay-triangle/remap path only if those new contracts are actually proved.
The objective remains t6-10 and the longest GIZMO duration.

[Dobrev, Kolev & Rieben](https://computing.llnl.gov/sites/default/files/SISC-12.pdf)
is background for compatible high-order kinematic/thermodynamic spaces;
its stability results are not imported into this mesh or remap.
