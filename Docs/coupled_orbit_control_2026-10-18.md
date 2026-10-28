# Gravity-only control qualifies coupled-response growth

Integrated the identical initial reference/plus/minus states under the
configured spline central gravity only, using independent RK4 at dt=.0005
and .00025 for duration .1. Preserved frozen exterior; verified no active
boundary membership change. No production trajectory or force was altered.

| Squared response final/initial | Hydro | Gravity only |
|---|---:|---:|
| Velocity | 1.289114 | 1.314078 |
| Circular-departure velocity | 1.082754 | 1.077569 |
| Omega-weighted displacement | 1.307331 | 1.327859 |

Gravity timestep central-velocity response relative difference is2.56e-9.
The observed increases are therefore also present without hydrodynamic
forces. Norm differences alone do not furnish additive causal fractions,
because the reference trajectories differ. The prior coupled response is
not evidence of a uniquely hydrodynamic growing mode.

Rebuilt Voronoi areas at the gravity-only endpoints and formed an
isentropic pressure proxy P=P_start*(V_start/V_end)^gamma for each case.
Initial areas agree within3.25e-14. On fixed initial r<1.9 IDs, acoustic
response ratio is1.543816 for this proxy versus1.496244 in actual hydro.
The proxy is not a hydrodynamic solution: pressure is never fed back into
the motion and no dissipative heating is applied. Nevertheless, growth
of the acoustic diagnostic alone also fails to establish self-amplification.

Interpretation: this seed mainly probes orbital/geometrical response that
already exists without fluid feedback. Do not interpret the31% combined
metric increase as a discovered instability or extrapolate it to t6-10.
Long-time physical failure remains real and unresolved.

Next priority: quantify continuously generated isentropic pressure forcing
along exact circular reference geometry versus the actual early t2/t4
geometry, using the same material IDs and linear pressure-force operator.
This targets forced departure directly instead of promoting an arbitrary
short-time norm to an instability diagnosis. Half-amplitude/hydro timestep
checks remain prerequisites if a later argument uses this particular
coupled tangent quantitatively; they are not yet completed.

Evidence: canonical_coupled_orbit_control_2026-10-18.json,
canonical_coupled_orbit_pressure_proxy_2026-10-18.json.
Scripts: audit_coupled_orbit_control.py, audit_coupled_orbit_pressure_proxy.py.
Remote root: /gpfs/kjhan/LagForce/Codex_coupled_orbit_control_20261018.
