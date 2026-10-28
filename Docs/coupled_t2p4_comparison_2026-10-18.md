# Continued response departs from gravity control by t2.4

Analyzed complete common checkpoint343 at t2.4004059177633943, duration
.3975 from seed184. Integrated gravity-only controls at this exact duration
with dt caps .0005/.00025; maximum squared-ratio change is4.65e-10.

| Squared response ratio from seed | Hydro | Gravity only | Hydro/gravity |
|---|---:|---:|---:|
| Raw velocity | 2.659966 | 2.621361 | 1.014727 |
| Circular-departure velocity | 1.846103 | .982941 | 1.878142 |
| Omega-weighted displacement | 2.798445 | 2.959017 | .945735 |

Hydro acoustic response ratio is3.707808; gravity control here has no
acoustic counterpart. Hydro even/odd velocity ratio is.001604. Relative
to gravity, the circular-departure squared response progressed from1.0048
at duration.1 to1.0842 at.1975 and1.8781 at.3975. This shows an emerging
difference between the fluid and gravity-only trajectories, obscured by
the raw velocity/position norms. It is not yet a converged instability
or an exponential eigenmode: hydro amplitude/dt tests remain running.

## Extension to t2.6004059

Complete checkpoint423 has duration .5975. Matching-time gravity comparison:

| Squared response ratio from seed | Hydro | Gravity only | Hydro/gravity |
|---|---:|---:|---:|
| Raw velocity | 4.305690 | 3.845076 | 1.119793 |
| Circular-departure velocity | 3.901254 | .962937 | 4.051414 |
| Omega-weighted displacement | 4.387236 | 4.307167 | 1.018590 |

Gravity dt comparison changes ratios by at most3.13e-10. Hydro even/odd
velocity ratio is.002512 and acoustic ratio6.314687. The departure response
is now materially larger than the gravity-only reference, while raw
position differences remain similar. This strengthens evidence of a fluid
response beyond orbital phase separation for this finite perturbation.
It still does not distinguish a transient response from an unstable mode,
or establish hydro timestep/amplitude convergence at this longer horizon.

Evidence: canonical_coupled_t2p6_comparison_2026-10-18.json. The half-amplitude
.1975 test has completed separately; fine-dt and restart controls remain
running. Refinement checks must eventually reach this longer interval;
do not transfer their shorter-horizon accuracy claims without measurement.

Both the t4 continuation(session62605) and independent .1975 convergence
batch(session46690) remain live at last check. Once the convergence batch
completes, verify it and continue its existing endpoints to this .3975
duration to test whether the larger difference persists under refinement.
Preserve all checkpoints; do not rerun the original branch from t0.

Evidence: canonical_coupled_t2p4_comparison_2026-10-18.json.
The historical t6-10 failure and long-term CUDA/OpenMP/MPI physics goals
remain unmet; these are early mechanism diagnostics.
