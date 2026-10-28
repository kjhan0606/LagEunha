# Exact-orbit kinematic volume control through t20

The canonical initial positions were evaluated at prescribed exact circular
orbits under the same softened gravity, r<2.5; exterior points were frozen.
Independent SciPy Voronoi geometry was constructed at t=0,2,...,20. Interior
r<1.9 cells are all finite and positive. Initial areas match saved geometric
volumes within 1.94e-14 relative. This is a geometry diagnostic: no simulated
trajectory was steered, no forces were integrated and no stability pass is
claimed. The result cannot replace the original hydrodynamic evolution.

With fixed initial cell masses and no entropy production, the pressure proxy
is P/P0=(V0/V)^gamma. Maximum sampled mass-weighted pressure-relative RMS is
.01588 for r0<.5, .006630 for .5<=r0<.6, and .002846 for .6<=r0<1.9. At t20,
inner volume ratio spans .94967 to 1.05256 and pressure-relative RMS=.01387.
The eleven samples do not bound extrema between samples or prove asymptotic
boundedness. Radial bins differ slightly from density-selected prior reports.

This qualifies the initial pressure-source finding: sizeable instantaneous
Vdot/V does not itself imply a large secular pressure increase. Along exact
orbits, geometric distortion produces much smaller sampled pressure errors
than the actual dissipative t12 trajectory. The open physical mechanism is
amplification of the initial seed through coupled pressure/velocity dynamics,
with entropy accumulating along the perturbed motion. Neither this control
nor the existing late-time QBP decomposition proves a unique feedback cause.

Next useful causal measurement is the response of the actual perturbation
equations to the measured geometric pressure seed, including pressure work,
gravity and dissipative damping. Do not replace the dynamics with prescribed
orbits or zero pressure-volume work to obtain an apparently stable result.

Evidence: canonical_circular_kinematic_volumes_2026-10-18.json.
Script: audit_circular_kinematic_volumes.py.
Remote output: /gpfs/kjhan/LagForce/Codex_circular_kinematic_volumes_20261018.
