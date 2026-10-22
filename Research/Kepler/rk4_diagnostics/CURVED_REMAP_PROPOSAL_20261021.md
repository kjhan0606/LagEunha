# Actual curved Kepler remap: proposal for one bundled Astra decision

This is a proposal, not an accepted algorithm or completed geometry gate.
Read MATERIAL_P3_RESULTS.md, MATERIAL_FINITE_RESULTS.md, the finite JSON and
the small native curved-map algebra contract before choosing the next bundle.

## Available exact structure

The actual current P3+two initial-velocity-residual map has degree<=6
power-plane numerators on fixed reference-radius rays. The denominator is
positive. Native tests verify positions, velocities, radial/angular
derivatives and J against the original map, including a truly curved state.
Use the actual packet, not chi=X+h*v0 alone or straight edge chords.

Proposed overlap: isolate all real roots of each active power-plane
polynomial within angular source-triangle intervals. Intersect the resulting
intervals with target halfspaces. Integrate original dm=rho0(r)*r dr dtheta
and transported actual p,K,L,physical U,area on the common refinement.
The radial reference density is exact, not reset current density. Physical
thermal energy is carried from P0*J^-gamma in this finite CONTROL; after
remap do not assign P0 or an initial analytic state to the target.

The screened weight relation remains M_i(W)+mu_i*W_i=m0_i. At a simple
active root the analytic mass derivative has nonnegative measure
rho0*r/(2*abs(dtheta[(qj-qi).chi])) dr and gives a symmetric face Laplacian.
Use actual roots and derivatives, not quantized ownership finite differences.

Hard issues requiring explicit scope/stop criteria: missed radial topology
events, tangencies/double roots, tiny closed lenses, uncertain polynomial
signs, global map injectivity, and boundary terms when active target changes.
An adaptive radial quadrature discrepancy alone cannot certify absence of a
tiny overlap. An edge-plane preimage is A3(s)+Omega(r(s))*B1(s)=0; fourth-power
elimination yields degree<=18 but adds roots that must be checked in the
original signed equation. This is an option, not an implemented solution.

The24 selected native cells comprise a bounded material patch with known
external traction, including disconnected pieces. Target power cells are
restricted to its actual curved fluid domain. This is not the full box.
At h=0,W=0 the original Voronoi partition is recovered. This identity alone
does not prove the user's dynamical W=0/individually-fixed-mass reduction;
its implementation remains an explicit requirement.

## Cold re-expression option to evaluate, not yet implement

Avoid compressing E then subtracting a changed K, and avoid an independent
kinetic reservoir. Seek actual conforming target velocity DOFs preserving
cell momentum, angular momentum and kinetic energy, while carrying physical
thermal energy/density separately and positively. Shared trace constraints
couple cells. Per-cell minimum K at given p,L is translation+rigid rotation,
but that local fact does not prove global conforming feasibility.

Possible construction: first find a conforming velocity satisfying linear
p,L constraints and per-cell K<=source K (convex quadratic inequalities).
Then use actual interior bubble directions with zero p/L and zero face trace
to restore each resolved kinetic energy. Minimize error to the transported
source velocity. If no admissible state or no resolved direction exists,
retain source fragments or refine; do not inject arbitrary noise, heat or a
separate reservoir. Global K conservation alone is insufficient for each
cell's actual energy/face ledger. Any claim of a face flux must be derived
from the overlap, not assigned afterwards to explain a residual.

Please decide whether this constrained projection is mathematically and
physically appropriate, or choose a simpler exact fragment representation
that can eventually run on the user's native CUDA/MPI path. Retaining
fragments is an intermediate representation, not an excuse for unbounded
memory or indefinitely avoiding an actual runtime. Preserve Re=infinity,
no steering, W0 reduction and the historical/longest-paper target.

## Production concerns

Current finite force/moment kernels already run on real CUDA and MPI, but
per-point host transfers and a root dense Schur remain reference costs.
Next production work needs device cell/triangle reduction, cached reference
mass/trace and sparse/matrix-free coupled solves. Keep physical algorithm
choice ahead of large optimization; no new time-integrator/degree sweep.
