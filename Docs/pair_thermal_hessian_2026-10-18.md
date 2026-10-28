# Negative thermal-energy curvature in a mixed pair displacement

Selected the nearest pair among generators at r<.6 in actual t10 geometry:
IDs 156/157, separation .00413165093, gap/mean sqrt(volume)=.16951.
At fixed per-cell mass and entropy, rebuilt independent Voronoi volumes for
center-of-mass-preserving relative displacement along pair-normal, tangent,
and their normalized sum. Included all r<2.4 cells in internal energy;
initial areas match production volumes within 3.42e-14 relative.

For U_i(V)=U_i0*(V_i0/V)^(gamma-1), directional curvature splits as

    U'' = sum gamma P/V * (V')^2 - sum P V''.

Normal curvature is +1.50676e-4, tangent +3.48649e-5, but mixed-direction
curvature is -1.84477e-4. Its compressibility term is +3.84688e-4 and
pressure/geometry term -5.69165e-4. The two-dimensional restricted Hessian
has a negative eigenvalue. Three perturbation sizes (.001,.0005,.00025
times gap) preserve this sign with small finite-difference variation.

Thus normal separation alone would miss a negative mixed deformation
direction. Nonuniform pressure and changing geometry can outweigh positive
compressibility stiffness in the evolved state. This is a local energy
landscape result, not proof of an exponentially growing full Kepler mode:
the state is moving and not stationary; gravity, rotation, damping, other
particle displacements and topology changes are not included. It does not
establish that this negative direction caused the original onset.

Next: verify pressure-force directional derivatives against the energy
curvature and compare constant-pressure/initial-entropy controls before
considering a geometry or pressure-force correction. Preserve conservation;
do not clip negative curvature or pressure work without a consistent model.

Evidence: canonical_pair_thermal_hessian_2026-10-18.json.
Script: audit_pair_thermal_curvature.py, default/--tangent/--diagonal.
Remote raw results: Codex_pair_thermal_curvature[_tangent/_diagonal]_20261018.
No existing trajectory modified.
