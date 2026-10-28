# Early-onset trajectory response at t2

Three independent restart branches completed from the preserved canonical
t2 checkpoint: reference and +/-0.001 times the naturally evolved velocity
departure from circular flow inside current r<0.6. Removed the seed's
mass-weighted translation. Only initial velocities were perturbed; no
steering was applied during evolution. Original long runs were preserved.

All branches exited zero with 40 accepted finite, zero-retry steps at
dt=0.0025, from t=2.0029059177634028 to t=2.1029059177634006.
Persistent particle IDs and masses match, and initial positions, volumes,
pressures and thermal energies are identical across branches.

| Diagnostic | Final / initial |
|---|---:|
| Reference inner circular-departure kinetic energy | 1.16059865 |
| Central velocity response norm | 0.97200633 |
| Central circular-departure velocity response norm | 0.96552796 |
| Central velocity-plus-acoustic quadratic diagnostic | 0.96174576 |

The final quadratic diagnostic consists of kinetic ratio 0.94479630 and
acoustic ratio 0.01694946. Even velocity response relative to odd response
is 4.53e-5 at this seed amplitude. This is a central finite perturbation
response, not an independently converged tangent solve. Half-amplitude and
half-timestep checks have not yet been performed.

The reference inner error grows while this additional velocity-only seed
damps. These measurements use different quantities and spatial supports:
reference energy uses initial r<0.5 and initial density<0.5; tangent norms
sum over all particles after a current-r<0.6 initial seed. Their numerical
ratios should not be directly subtracted or treated as an energy budget.

This limits the hypothesis that the observed early growth is immediate
self-amplification of this velocity pattern alone. It does not rule out
collective growing modes involving position, density or entropy, transient
growth on longer intervals, or ongoing forcing by geometric pressure and
viscous truncation errors. The positive velocity/acoustic diagnostic is not
a conserved perturbation energy and omits displacement and entropy degrees.

Next priority: attribute the growing reference departure over this same
interval to pressure, differential rotation, TCF and split HG using accepted
trajectory data, rather than inferring forcing from the damped tangent alone.
Then test the implicated coupled state perturbation. Do not adopt a pressure
or heat suppression fix, or claim long-term physical stability from this test.

Evidence: canonical_onset_velocity_tangent_2026-10-18.json.
Runner: run_onset_velocity_tangent.py.
Analyzer: analyze_onset_velocity_tangent.py.
Remote root: /gpfs/kjhan/LagForce/Codex_onset_velocity_tangent_20261018.
