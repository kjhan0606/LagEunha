# Continuum energy bound versus the native grid-scale shear disturbance

Purpose: separate physical KH amplification from numerical amplification
without fitting a growth rate, adding viscosity, or launching another long
trajectory. Uses saved zero-seed32x64 and64x128 states only. This is an
analytical/postprocessing diagnostic, not a replacement Python simulation.

## Derivation for the actual smooth initial profile

Background is rho0(y), velocity (U(y),0), constant pressure P0=10 and gamma.
Write Eulerian perturbations u,v,p. The closed linear velocity/pressure system
is

    (partial_t + U partial_x) u = -U' v - partial_x p/rho0,
    (partial_t + U partial_x) v = -partial_y p/rho0,
    (partial_t + U partial_x) p = -gamma P0 (partial_x u+partial_y v).

Density/entropy perturbations do not enter these three equations because
background pressure is constant. On the periodic domain define

    E = (1/2) integral [rho0 (u^2+v^2) + p^2/(gamma P0)] dV.

Pressure work cancels by integration by parts, background x-advection has
zero integrated work, and

    dE/dt = -integral rho0 U' u v dV,
    |dE/dt| <= max|U'| E.

For U=-1+tanh((y-.5)/.05)-tanh((y-1.5)/.05), max|U'| <=20.
Therefore E(t2)/E(t1) <= exp[20(t2-t1)] for the unforced linear system.
This bounds arbitrary superpositions, not just eigenmodes or a single k.
It is an ENERGY bound; the corresponding norm/amplitude exponent is10.

Ten independently generated periodic smooth test fields verify the energy
identity using skew-adjoint Fourier derivatives, max absolute residual
6.67e-16. This checks the factors/signs, not native solver consistency.

## Applying the diagnostic to saved native states

Use u=vx-U(y), v=vy, p=(gamma-1)Uthermal/V-10. Integrate with actual cell
volumes and rho0 evaluated at generators. This is center quadrature on a
moving mesh, not a rigorous continuous reconstruction. A second norm replaces
rho0 V by the stored mass; both estimates are close (well below1% difference
in the small-amplitude states), too small to explain the discrepancies below.

| Nx | interval | measured E ratio | continuum upper ratio | ratio / upper bound |
|---|---|---:|---:|---:|
|32|.6--1|3.86677e6|2.98096e3|1297|
|32|1--1.5|1.39700e8|2.20265e4|6342|
|64|.6--1|2.52786e9|2.98096e3|848002|

All three endpoints have max velocity perturbation <6.43e-4 and max relative
pressure perturbation <1.92e-4. They precede the large nonlinear64,t=1.5
disturbance; that later interval is explicitly excluded from a linear bound
claim. Observed logarithmic energy growth rates are approximately37.9,37.5,
54.1, versus the bound20.

At very early times numerical errors may continuously seed perturbations,
so exceeding the homogeneous bound need not mean every perturbation follows
a single unstable discrete eigenmode. Either excess amplification or ongoing
forcing is numerical content, not unforced continuum KH growth. Center
quadrature and finite linearization remain qualifications. Taken with the
two-cell dominant wavelength and timestep-refinement control, the discrepancy
is strong evidence that the native response cannot be qualified as physical
KH growth merely because its mass/entropy are conserved.

## Development implication

The existing successful Kepler force/volume-work path is useful but insufficient
for a common hydro method. Do not hide this result with a smoothing/viscosity
switch. The next discriminating spatial calculation is the linear pressure
response along the moving shear geometry: uniform-pressure force cancels,
but its volume Jacobian changes as rows pass one another. Compare its
short-wave response with this continuum bound before choosing a correction.
This audit alone does not prove a particular missing term or justify a new
physical variable. Existing1D reference, dynamic-weight work, and2D trajectories
remain preserved;3D is still deferred and GIZMO runtime qualification is open.

Artifact: Research/MaterialVolume_20261105/kh_continuum_growth_bound_20261106_v1/
audit.py and results.json. No native source or physical parameter changed.
