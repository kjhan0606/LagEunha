# Current KH comparison with the existing linear reference — 2026-11-09

This is a diagnostic audit of completed 2D runs, not a new hydro method or an accuracy pass. No production source, live run, dissipation coefficient, or initial condition was changed. All twelve native saved states were re-extracted from checkpoints and matched byte-for-byte after sorting by ID to their archived state arrays. Checkpoint times and input hashes are recorded in the reproducible report.

## Early growth

The reference is the existing independent compiled-C linear Euler calculation (`kh_linear_reference_20261106_v1/linear.c`), with the same continuous smooth KH background and imposed x fundamental. Native total transverse kinetic energy includes all Fourier modes; reference kinetic energy includes only the linear imposed mode.

| Run | t=0.3 relative difference | t=0.6 relative difference |
|---|---:|---:|
| Current 128×256, C=0.2 | −0.8282% | −2.0762% |
| Previous 128×256, C=0.2 | −0.7966% | −2.1143% |
| Current 64×128, C=0.8 | −3.1996% | −10.6153% |
| Previous 64×128, C=0.8 | −2.9307% | −10.1632% |

The 512→1024 reference grid change is below 5.4e−8 relative at these times; halving its timestep changes it by less than 1.9e−11. Reference numerical precision is therefore much smaller than these diagnostic discrepancies. Finite-volume sampling and nonlinear/higher-mode contributions remain interpretation limits. Different CFL coefficients in the two native resolutions prevent calling this a pure spatial-convergence result.

## Late fundamental, t=1.35

The existing spectral linear extrapolation gives amplitude 0.0970358266893. Native current 128×256 gives 0.0906793433600; previous 128×256 gives 0.0906831479167. Their difference is -0.004195 percent relative to the previous amplitude. Synchronized GIZMO128 gives 0.0808276348322. Current native64 gives 0.0775468373250 and synchronized GIZMO64 gives 0.0500345929944.

These are point-mass Fourier integrals of transverse momentum divided by total mass, not smoothed density images. Late linear extrapolation is not an exact finite-amplitude nonlinear solution, so closeness to it does not establish that either hydro method is more accurate. Both GIZMO comparisons use synchronized endpoint velocities; the existing input audit still reports approximately 0.058% pressure offset after GIZMO initialization.

## Consequence for development

Matching the endpoint volume and pressure work did not substantially change the resolved KH fundamental. Thus the remaining native–GIZMO growth discrepancy cannot be attributed primarily to that recent correction on this evidence. Keep the geometric consistency fix; do not tune another viscosity/sensor to fit one late linear curve. Complete the already-running fixed-source Sedov refinement and RT evolution, and retain the matched-state timestep controls separately from resolution/physical-accuracy claims. The ≤2× GIZMO requirement remains unproven at matched accuracy.

## Artifacts

- `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_kh_linear_20261109_v1/audit.py`
- Same directory: `comparison.json` (hashes, early/late values, reference refinement, limitations), `growth.png` (inspected), and re-extracted states.

No 3D evolution or manuscript claim is introduced.
