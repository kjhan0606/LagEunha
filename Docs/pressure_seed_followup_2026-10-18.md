# Distributed pressure seed: duration and timestep follow-up

Resumed existing +/- pressure branches from checkpoint 2313 for 32 further
steps, giving duration .020. Reused the completed same-binary reference at
the identical endpoint. Separately evolved +/- seeds and unperturbed reference
from the original checkpoint with half dt=.00025 for 16 steps to duration
.004. All five runs exited zero; all expected accepted rows are finite and
zero-retry. No original long trajectories were restarted or changed.

| Duration | dt | Kinetic ratio | Acoustic ratio | Combined ratio |
|---|---:|---:|---:|---:|
| .004 | .0005 | .001424405 | .999801417 | 1.001225822 |
| .004 | .00025 | .001424521 | .999801797 | 1.001226318 |
| .020 | .0005 | .026704008 | .972715505 | .999419512 |

Ratios normalize the positive velocity/pressure quadratic diagnostic to its
initial pressure-seed value. Short-interval normalized response difference
between timesteps is 1.59e-5; combined ratios differ by 4.95e-7. The small
short-interval rise is reproduced with half dt but does not persist at .020.
This seed therefore shows pressure-to-velocity response without established
sustained amplification over the sampled intervals. It is not evidence of a
stable full solver, nor does it identify repeated geometric forcing as the
unique remaining cause. Longer-interval timestep convergence was not tested.

Together with the damped local velocity seed, these results limit the use of
the frozen two-particle negative-curvature model as a diagnosis of actual
growth. A causal instability claim now needs an appropriate collective or
position/entropy perturbation, or direct attribution of continuously generated
errors along the existing trajectory. Do not declare success from near-unity
short diagnostic norms while historical t6-10 physical failure persists.

Evidence: canonical_pressure_seed_followup_2026-10-18.json.
Scripts: run_pressure_seed_followup.py, analyze_pressure_seed_followup.py.
Remote outputs: /gpfs/kjhan/LagForce/Codex_pressure_followup_20261018.
