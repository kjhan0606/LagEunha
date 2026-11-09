# Gresho radial-error localization

Saved native 32x32 checkpoints, frozen and dynamic weights, t=.015,.03,.06,.1.
Same dt cap 1e-4 and restart endpoints; no evolution rerun or physics change.
Artifact: Research/MaterialVolume_20261105/gresho_radial_audit_20261106_v1.
extract.c uses the native checkpoint codec. analyze.py checks 1024 finite
rows per state. radial_metrics.json stores volume-weighted radial diagnostics.

At t=.1 the dynamic-weight radial RMS is .010168975 versus .006902235 frozen.
The annulus .18<=r<.22 contains 85.10% of dynamic radial squared error;
frozen contains 26.86%. Its conditional RMS is .04528081 dynamic versus
.01726066 frozen. Every other tabulated radial interval has a smaller
conditional RMS with dynamic weights. Thus the global worsening is localized
around the r=.2 velocity-gradient break, not uniformly throughout the shear
or rigidly rotating core.

Dynamic annular mean radial velocity is +.00628675, +.01331821, +.02672392,
+.04221361 at the four times. The outward drift is growing over this saved
interval; it is not merely a larger zero-mean scatter in the endpoint plot.
Frozen annular means are +.00424222,+.00241809,-.00530299,-.01348475.
These are spatial annular selections at each time, not tracked material IDs.

This is evidence against adopting dynamic weights based only on density L1
or energy conservation. It does not establish the cause: pressure-force
truncation, weight dynamics and point-sampled initial balance must still be
separated. The checkpoints store accelerations, but their stage/freshness
has not been established, so no causal acceleration attribution is made.
Next targeted comparison should evaluate the native RHS on the same saved
states with geometry and thermodynamics held fixed, before adding damping
or another empirical term. Keep frozen successful-Kepler ingredients as the
control; retain all dynamic results. Full 2D KH/RT, shock compatibility with
1D, and GIZMO end-to-end timing remain incomplete. Defer 3D.

## Native RHS replay and controlled rate-zero comparison

All eight saved states were passed through the current native C/CUDA K1
probe on Syntax A10 (one process, two OpenMP threads). Each exited zero
with the computed-K1 marker; no accepted evolution. Source checksums and
environments are saved per rhs_MODE_TIME directory. The diagnostic snapshots
contain freshly written GPU accelerations, not stale endpoint accelerations.

For dynamic states, annular mean a_r + v_phi^2/r is +.46540,+.46601,+.43179,
+.31716 at t=.015,.03,.06,.1. This is consistent with the measured outward
radial drift, though the spatial annulus does not track identical particles.

At t=.1 a counterfactual diagnostic zeros canonical weight momentum pi,
leaving positions, velocities, W, mass, U and V unchanged. The native K1
acceleration is exactly identical (maximum component difference zero).
Thus instantaneous Wdot is not directly producing this pressure-force error.
This does not exclude the accumulated effects of earlier weight evolution.

A second diagnostic also replaces U by P_exact(x)*V/(gamma-1). The annular
outward acceleration residual falls from .31716181 to .19598887, but remains
positive. This is point-sampled analytic pressure, NOT the exact cell-average
pressure used by the initial quadrature. It is not a proposed energy reset,
not a valid conservative evolution step, and not proof that geometry alone
causes the remainder. No altered snapshot is an accepted restart.

Decision: do not add a Wdot limiter to fix the measured radial force; this
same-state test supplies no direct force mechanism for that intervention.
Next useful comparison is a consistent cell-integrated equilibrium/pressure
force check near the velocity-gradient break, retaining the native W0
control. Pressure averaging, force truncation and accumulated geometry
changes must be distinguished before another evolution-model modification.
Artifacts: run_rhs.py, rhs_results.json, analyze_rhs.py, rhs_metrics.json,
counterfactual.c, run_counterfactual.py, counterfactual_results.json.

## Cell-integrated analytic pressure diagnostic

Used the exported native power-cell faces at t=.1 for both 32x32 states.
Triangular Duffy quadrature of orders 12,24,48 integrates analytic Gresho P
within each actual cell. Maximum pressure-average change 24->48 is
7.89e-9 frozen and 1.03e-8 dynamic. Summed triangle area agrees with native
V within 3.26e-19. A single numerically zero-length dynamic face produced
signed area -2.71e-20; skipped only negative cross products >=-1e-18,
with the full-cell area agreement explicitly checked. Nondegenerate negative
triangles are rejected, not silently repaired.

Diagnostic snapshots set pi=0 and U=P_cell_average*V/(gamma-1), preserving
positions, mass, W and velocity. Native K1 completes for both; no evolution.
These are forced-pressure counterfactuals, not proposed thermodynamic resets.
They intentionally break the trajectory's entropy history to isolate the
force response. Quadrature is offline and is not a production algorithm.

Annular means .18<=r<.22:

| State | original a_r+v_phi^2/r | cell-average analytic P replacement |
|---|---:|---:|
| frozen | -.26779891 | +.33539891 |
| dynamic | +.31716181 | +.23138050 |

Using exact Gresho v_phi at the same points instead of evolved velocity gives
+.40753559 frozen and +.27413448 dynamic for the pressure-replaced states.
Thus pressure averaging alone does not restore the point-orbit equilibrium
in either geometry. The dynamic path is not uniquely responsible for this
residual. This narrows the investigation to the spatial pressure-force and
representation of orbital balance, with evolved density and positions still
held fixed; it does NOT isolate a geometry-only error or prove a defect in
the variational derivation. In particular, rho is not reset to unity.

Before altering the model, compare the initial uniform-density geometry
with this saved-state result and examine the spatial pressure-force order.
Do not cure the error by assigning analytic pressure during evolution, adding
an empirical radial force, or selecting the method using density L1 alone.
Artifacts: integrate_pressure.py, quadrature_results.json, cell_pressure.c,
run_cell_pressure.py, cell_pressure_results.json. All probes terminal rc0.
