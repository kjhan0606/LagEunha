# Actual first-step pressure response

Two isolated branches resumed the preserved t0 checkpoint with identical
initial positions, velocities, pressure, volume and internal energy to the
canonical initial snapshot. Both executed one accepted HG/RK4/HG step,
dt=1e-4 or 5e-5, with no retries and finite state. Original trajectories were
not restarted or changed. These branches test an initial derivative only.

Pressure secants (P_end-P_start)/dt were compared against the independently
saved t0 expression -gamma P Vdot/V+(gamma-1)(H_TCF+H_HG)/V.

| Region (initial snapshot radius) | Relative L2 error dt=1e-4 | dt=5e-5 |
|---|---:|---:|
| r<.5 | .00259463 | .00129732 |
| .5<=r<.6 | .00275854 | .00137927 |
| .6<=r<1.9 | .00205231 | .00102616 |

The errors halve with dt, as expected for a forward secant approaching an
initial derivative. This is NOT an RK4 convergence-order test. Inner measured
relative pressure-rate RMS at the smaller step is .1985948 versus predicted
.1986028. Heating alone gives approximately unit relative error; omitting
compression fails to explain the initial pressure change. Compression alone
has about .01190 relative error in that region, consistent with a smaller
but nonzero dissipative contribution.

These findings validate the initial pressure-seeding pathway in actual
evolution, beyond an algebraic volume-rate observation. They do not show that
this pathway alone causes t6-10 collapse or justify setting volume derivatives
to zero. Next physical requirement is a conservative spatial correction that
reduces this equilibrium defect while respecting no steering and the approved
weight/mass reduction; it must then survive the historical failure interval.

Initial-region masks here use the checkpoint coordinates; prior instantaneous
tables use reconstructed RHS coordinates, so boundary-radius cell membership
can differ at roundoff near r=.5. Do not compare their regional RMS as if the
selected IDs were necessarily identical.

Evidence: canonical_initial_pressure_step_2026-10-18.json.
Analyzer: analyze_initial_pressure_step.py.
Remote runs: /gpfs/kjhan/LagForce/Codex_initial_pressure_step_20261018.
