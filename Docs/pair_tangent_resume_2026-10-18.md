# Resumed pair response through total duration .020

Continued the four existing +/- velocity-seed branches and same-binary
reference from their final checkpoint 2313. No seed was reapplied and no
original long trajectory was restarted. Five runs completed 32 more actual
steps, giving 40 steps total from t10.0020667693911 to 10.022066769391124.
All accepted rows are finite with zero retries.

The central velocity-response mass norm relative to its initial value is
.484348 at duration .020, compared with .904062 at duration .004. Half/full
seed amplitudes agree to 5.11e-7 in the full velocity derivative. About 4.17%
of squared mass-weighted velocity response is now outside pair 156/157.

Also measured the positive quadratic diagnostic

    sum m |delta v|^2 + sum V (delta P)^2/(gamma P)

using unperturbed endpoint P,V. Divided by initial sum m|delta v|^2, it is
.250659 at .020 (acoustic part .0160660); at .004 it was .818684. This norm
does not include displacement/tidal terms or entropy degrees of freedom and
is not a conserved perturbation energy for the evolving shearing flow.

For this seed and interval, the actual solver strongly damps velocity and
kinetic/acoustic response despite the restricted pressure-only negative
curvature and frozen-pair positive growth eigenvalue. The restricted model
is not a reliable quantitative prediction of this actual response. This does
not exclude a position/entropy seed, collective mode, nonnormal transient
growth or later growth. It is not a long Kepler stability pass.

Next physical priority should move from this damped local velocity seed to
spatially distributed pressure/volume excitation and its repeated forcing.
Preserve all branch endpoints for any justified continuation.

Evidence: canonical_pair_tangent_resume_2026-10-18.json and
canonical_pair_tangent_acoustic_2026-10-18.json.
Scripts: resume_pair_tangent_branches.py, analyze_pair_tangent_resume.py.
Remote outputs: /gpfs/kjhan/LagForce/Codex_pair_tangent_resume_20261018.
