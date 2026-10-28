# Initial geometric compression seeds pressure structure

Using existing t0 RHS, common-face geometry and instantaneous reconstructed
HG heat, fixed-mass ideal-gas pressure obeys

    Pdot = -gamma P Vdot/V + (gamma-1)(H_TCF+H_HG)/V.

Initial pressure for r<1.9 is uniform to floating precision, about
2.777779632144756e-5. The independent pre-HG energy identity
(gamma-1)(Edot-v.F)/V-P Vdot/V agrees with compression+TCF pressure derivative
to 9.17e-14 absolute in that region, limited by saved geometry precision.

Mass-weighted RMS of relative pressure derivatives:

| Current radius | compression | TCF heat | HG heat |
|---|---:|---:|---:|
| r<.5 | .201479 | .00259912 | .0000120608 |
| .5<=r<.6 | .0864307 | .0000702418 | .0000119271 |
| .6<=r<1.9 | .0220130 | .00000692730 | .000000619968 |

Applying the linear pressure-force operator to each Pdot at fixed geometry
gives inner pressure-acceleration derivative RMS .0170833 from compression,
8.16e-5 from TCF heating, and 9.34e-7 from HG heating. These are pressure-only
derivatives, not total jerk; direct dissipative forces, gravity, subsequent
geometry changes and nonlinear feedback are not included. Small initial
constant-pressure residual in reconstructed geometry is retained, not hidden.

This provides a quantified initial pathway: correct geometric volume change
under circular shear produces density/pressure structure even though material
divergence is zero. At t0 the compression contribution dominates pressure
seeding; at t10 evolved velocity perturbations dominate QBP signals and TCF
dominates sampled entropy production. These two observations are compatible.
They do not yet prove this pathway accounts for the entire long-time failure.

No trajectory or initial condition was changed. No artificial pressure or
volume suppression is proposed: any correction must preserve pressure-work
compatibility, conservation, no steering and the approved weight reduction.
Next causal test should compare this initial derivative with sufficiently
early saved evolution and quantify the pressure-induced acceleration growth,
before adopting a modified spatial operator.

Evidence: canonical_initial_pressure_seed_2026-10-18.json.
Script: audit_initial_pressure_seed.py.
Remote data: /gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018.
