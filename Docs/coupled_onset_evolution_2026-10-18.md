# Coupled isentropic position/velocity response at onset

Evolved all three prepared branches with MPI1, OMP4 and dt=.0025 for40
steps, t=2.0029059177634028 to2.1029059177634006. All exited0 with finite
accepted states and no retries. Particle IDs/masses match at both endpoints.
No new steering or entropy adjustment was applied during evolution.

Central response uses +/-0.001 seeds from actual early position/velocity
departure, with isentropic thermodynamics reconstructed at the initial
perturbed geometry. Final/initial squared diagnostic ratios:

| Quantity | Ratio |
|---|---:|
| Velocity response, sum m*delta_v^2 | 1.289114 |
| Circular-departure velocity response | 1.082754 |
| Acoustic pressure response | 1.496244 |
| Omega-weighted displacement response | 1.307331 |
| Departure velocity + acoustic + displacement | 1.311700 |

Even velocity response is0.001066 of the odd response. The reference
global departure-energy ratio is1.150196, a different observable from
the central perturbation response. Prior velocity-only seeds damped;
this coupled initial perturbation grows over the measured interval.

This is evidence that the velocity-only probe omitted a relevant coupling.
It is not yet evidence of an exponentially unstable eigenmode. Displacement
and raw velocity responses can grow through orbital shear/phase separation;
positive norms are not conserved perturbation energies, and entropy was
not independently perturbed. No half-amplitude/time-step check has yet
been completed. The acoustic and circular-departure responses must be
retained separately when comparing controls.

Next: half-amplitude central branches using identical prepared geometry
rules, then timestep refinement and a noninteracting circular-orbit
variational control to quantify kinematic phase growth. Reuse the completed
reference where identical. Do not change the physical equations based on
one short response. Historical t6-10 failure and long evolution remain gates.

Evidence: canonical_coupled_onset_evolution_2026-10-18.json.
Runner: run_coupled_onset_evolution.py.
Analyzer: analyze_coupled_onset_evolution.py.
Remote root: /gpfs/kjhan/LagForce/Codex_coupled_onset_seed_20261018.
