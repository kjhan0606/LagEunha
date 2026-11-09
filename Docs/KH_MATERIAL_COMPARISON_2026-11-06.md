# New material KH versus prior fixed-weight baseline

New native constrained-weight, reconstructed acoustic work/force, SSPRK2 trajectory:32x64, no seed; compared with archived native W=0, fixed-mass, variational-pressure KDK at SAME32x64 and time0.6. Both have P=10, same smooth KH density/velocity profiles. This is a comparison of algorithm packages, not an isolation of weight effects or integrator effects.

At0.6 prior unseeded transverse kinetic energy1.23457e-23, relative pressure RMS6.32071e-13; new transverse kinetic3.17909e-7, relative pressure RMS2.10187e-4. Therefore the new method is worse at preserving this stationary shear in this interval, despite total-energy conservation. Do not describe conservation or reaching target time as an accuracy pass. Prior method later develops grid-scale modes, especially64x128; that does not erase its superior early stationary-shear result at32x64.

New row-decomposition at0.02,0.1,0.3,0.6 verifies exactly32 particles in each y-bin. The mass-weighted within-row transverse kinetic fraction is at most3.1e-19; essentially all new transverse motion is row mean (x-independent). This differs from the old later nonzero-x high-k mode. The direct Fourier projection is not exact Parseval on a deformed mesh, and the tiny nonzero-mode dominant k is not meaningful when its energy is negligible. Pressure/heat-driven row response is a hypothesis supported by accumulated acoustic heat; causal term isolation remains outstanding.

Artifacts:
- Research/MaterialVolume_20261105/material_kh_pilot_20261106_v1/history_metrics.json
- Research/MaterialVolume_20261105/material_kh_pilot_20261106_v1/unseeded_modes.py and unseeded_modes.json
- Research/MaterialVolume_20261105/kh_unseeded_20261106_v1/metrics.json

Continue existing accepted trajectories through1.5 to compare the previously problematic interval. Do not restart them. No new viscosity switch/physical parameter added. Low-pressure shear, RT,3D and GIZMO runtime target remain unqualified.
