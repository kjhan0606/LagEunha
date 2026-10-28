# Independent coupled response convergence checks: running

Target is the completed coarse checkpoint263 time2.2004059177633986,
duration .1975 from seed184. Prepared half-amplitude states with the same
isentropic geometry procedure; max area change1.875e-5 and entropy relative
error4.45e-16. All masses unchanged and reference source byte-identical.

Launched five independent MPI1 OMP4 cases with the same orbit-work binary:

| Case | Epsilon | dt cap | Expected accepted steps | PID |
|---|---:|---:|---:|---:|
| amplitude_plus | +.0005 | .0025 | 79 | 3114366 |
| amplitude_minus | -.0005 | .0025 | 79 | 3114373 |
| fine_reference | 0 | .00125 | 158 | 3114380 |
| fine_plus | +.001 | .00125 | 158 | 3114387 |
| fine_minus | -.001 | .00125 | 158 | 3114394 |

Root: /gpfs/kjhan/LagForce/Codex_coupled_convergence_20261018.
Runner session46690. Original t4 continuation remains session62605.
Do not restart a case on observation timeout. Completion requires actual
exit codes, final snapshots and analysis, not this launch record.

Prepared analyze_coupled_onset_convergence.py to check accepted steps,
finite/zero-retry states, times, IDs/masses, then compare full central
position/velocity/circular-departure/pressure response fields in weighted
L2. It also records common-final-weight growth ratios, which differ from
the earlier ratios with time-dependent weights.

The completed coarse branch was restarted after .1; new cases are
uninterrupted. Differences can include restart effects. If the discrepancy
is material, a matching coarse uninterrupted control is needed before
attributing it entirely to amplitude or timestep. No convergence or
long-time stability pass has yet been established.

Scripts: prepare_coupled_onset_seed.py now accepts --root/--epsilon;
run_coupled_onset_convergence.py; analyze_coupled_onset_convergence.py.

## Half-amplitude cases completed

Both amplitude cases exited0 with79 finite, zero-retry steps. Ran the
analyzer with --amplitude-only while fine-dt cases remain running.

| Normalized response field | Relative L2 difference | Squared norm ratio |
|---|---:|---:|
| Position | .000453823 | 1.00005230 |
| Velocity | .000848915 | 1.00008501 |
| Circular-departure velocity | .00751983 | 1.00056409 |
| Pressure | .00167513 | 1.00006822 |

Norms are close, but circular-departure field shapes differ by .752%.
This is not a bitwise pass or an asymptotic convergence proof. Restart
segmentation differs between cases, so launched matched full-amplitude
uninterrupted plus/minus cases at dt=.0025 to isolate that effect before
interpreting the amplitude difference. No completed cases were overwritten.

Restart-control root: Codex_coupled_restart_control_20261018.
Session28065; runner run_coupled_restart_control.py; analyzer
analyze_coupled_restart_control.py. Fine-dt cases remain session46690;
t4 continuation remains62605. After these checks, continue the existing
refined checkpoints to the already observed t2.4004059 comparison time.
Evidence: canonical_coupled_half_amplitude_2026-10-18.json.

## Half-timestep cases completed; refined continuation started

All three fine-dt cases exited0 with158 accepted finite, zero-retry steps.
At t2.2004059, half-dt central response relative L2 differences versus
coarse are .105% position, .196% velocity,1.552% circular-departure velocity,
and .548% pressure. Squared circular-departure norm ratio is .99897568
(about .1024% lower). Thus scalar norms are close while spatial response
has percent-level sensitivity. Do not claim formal convergence order or
spatial-field convergence from these two settings.

Evidence: canonical_coupled_convergence_2026-10-18.json. Restart control
session28065 is still pending at last check and must be used to isolate
its contribution to these differences.

Continued all five completed refined cases from their existing endpoints
t2.2004059 to t2.6004059, retaining .2 output intervals. No original
trajectory was restarted from its beginning. New root:
/gpfs/kjhan/LagForce/Codex_coupled_refined_to_t2p6_20261018.
Session93108. PIDs amplitude_plus3129961, amplitude_minus3129968,
fine_reference3129975, fine_plus3129982, fine_minus3129989. Coarse-amplitude
cases require160 further steps; fine-dt cases320. Runner:
resume_coupled_refined_to_t2p6.py. After actual completion analyze using:

    analyze_coupled_onset_convergence.py --root /gpfs/kjhan/LagForce/Codex_coupled_refined_to_t2p6_20261018 --coarse-step 423

The analyzer now accepts a root/coarse-step and source restart indices.
It still requires completed exit codes, correct step counts and matching
target times. These launch details are not completion evidence.

## Restart confound resolved for the completed central-response comparison

Session28065 completed both full-amplitude uninterrupted cases, exit0,
79 finite zero-retry steps. Analyzer verifies at t2.2004059 that all reader
fields id,x,y,vx,vy,mass,volume,den,pressure,ie,w2 are bitwise equal to
the40+39-step restarted plus/minus cases. Therefore the measured central
response differences at this endpoint are not caused by that restart.
This is not a claim about every snapshot field, unperturbed reference
restart, other settings, or later restart horizons.

Evidence: canonical_coupled_restart_control_2026-10-18.json.
The refined extension session93108 was verified live, with amplitude_plus
accepted step21 (elapsed .0525 from t2.2004059), finite and zero retries.
Original t4 continuation remains session62605; neither is complete yet.

## Refined intermediate checkpoint at t2.4004059

Both half-amplitude branches saved checkpoint343 and continued running.
Added --intermediate analysis mode: requires the saved time and the
corresponding accepted finite, zero-retry prefix, but does not require a
terminal exit code. Output filenames identify the analyzed coarse checkpoint.
This avoids labelling a live run complete or waiting unnecessarily for
the full .6 target.

Compared half-amplitude central responses at duration .3975:

| Field | Relative L2 difference | Squared norm ratio to full amplitude |
|---|---:|---:|
| Position | .00100177 | 1.00014099 |
| Velocity | .00141938 | 1.00009920 |
| Circular-departure velocity | .01090029 | 1.00035358 |
| Pressure | .00201889 | 1.00009783 |

The departure squared growth is1.846756 versus1.846103. Hence the observed
growth persists at half amplitude, although field shapes differ at roughly
one percent. This is evidence of amplitude robustness of the norm, not a
formal derivative convergence rate. Fine-dt branches had not yet reached
this checkpoint when the amplitude analysis ran. Do not claim temporal
convergence for this horizon until their actual outputs are compared.
Later restart histories differ; the earlier bitwise restart result is
limited to its tested t2.2004059 endpoint.

Evidence: canonical_coupled_t2p4_half_amplitude_2026-10-18.json.
Command:

    analyze_coupled_onset_convergence.py --root /gpfs/kjhan/LagForce/Codex_coupled_refined_to_t2p6_20261018 --coarse-step 343 --intermediate --amplitude-only

## Same t2.4004059 checkpoint: fine timestep now available

All three fine branches saved checkpoint502 at the same physical time.
The intermediate analyzer verified the160-step accepted finite/no-retry
prefix since their restart; runs continue toward t2.6004059.

| Half-dt response field | Relative L2 difference | Squared norm ratio to coarse |
|---|---:|---:|
| Position | .00229698 | .99967811 |
| Velocity | .00268116 | .99981068 |
| Circular-departure velocity | .01709378 | 1.00050225 |
| Pressure | .00634821 | .99917545 |

Departure squared growth is1.847030 for half dt versus1.846103 coarse and
1.846756 half amplitude. Growth of this scalar diagnostic persists under
both refinements with roughly .05% norm differences at this horizon.
Percent-level spatial response differences remain; this is not a formal
order estimate, a full mode-convergence claim, or a long-time stability pass.
The later restart-history caveat remains scoped as previously documented.

Evidence: canonical_coupled_t2p4_convergence_2026-10-18.json.
Same intermediate command without --amplitude-only analyzes both checks.
Maintain sessions93108 and62605, and next compare the completed t2.6004059
outputs using coarse checkpoint423.

## Completed half-amplitude check at t2.6004059

Both amplitude branches exited0 after160 accepted finite steps without
retries. At the matching checkpoint423, the circular-departure squared
growth is3.90125436 coarse versus3.89853860 half amplitude. The squared
norm differs by -0.06961%; the response field relative L2 difference is
1.42177%. Position, velocity, and pressure relative L2 differences are
0.18128%,0.22489%,0.31386%, respectively. This supports persistence of
finite-time growth under amplitude reduction, not an eigenmode proof.
The fine-dt branches are still running; their final convergence is pending.
Evidence: canonical_coupled_t2p6_amplitude_2026-10-18.json.

The unchanged coarse continuation also has complete common snapshots at
t3.0004059,t3.2004059,t3.4004059. Departure squared growth from the t2 seed
is16.4596,30.6606,79.6249, respectively. At t3.4 the even/odd velocity ratio
is0.007855. These later endpoints have NOT passed amplitude/dt comparison
or matching-time gravity control yet. They indicate worsening response,
not long-time stability. Preserve the ongoing continuation to t4.
Evidence: canonical_coupled_t3p4_response_2026-10-18.json.

## Completed full comparison at t2.6004059

All five refined branches exited0. Final analyzer verified160 accepted
steps per amplitude branch and320 per fine-dt branch, finite with no retries.
Departure squared growth: coarse3.90125436, half amplitude3.89853860,
half dt3.89364349. Half-dt squared norm is0.99804912 times coarse (-0.1951%).
Half-dt field relative L2 differences: position0.36491%, velocity0.36842%,
circular departure1.82879%, pressure0.73865%. Scalar growth persists, while
spatial percent-level discrepancies and restart-history caveats remain.
Evidence: canonical_coupled_t2p6_convergence_2026-10-18.json.

Continue from these completed endpoints, not the original seed, to the
existing coarse checkpoint743 time3.400405917763373. New root:
/gpfs/kjhan/LagForce/Codex_coupled_refined_to_t3p4_20261018.
Runner: resume_coupled_refined_to_t2p6.py with --source the completed t2p6
root, --root the new t3p4 root, --coarse-step743. Session23450.
Amplitude sources423 need320 steps; fine sources662 need640 steps.
Intermediate comparisons use the existing analyzer with --intermediate
and coarse checkpoints503,583,663; final uses743 without --intermediate.
Preserve original t4 continuation62605. These remain mechanism convergence
tests, not replacements for the historical long-time failure gates.

## t2.8004059 intermediate comparison completed

All five continuing refinement branches saved matching checkpoints:
amplitude503 and fine822. Accepted prefixes are finite with no retries.
Departure squared growth: coarse8.3526721, half amplitude8.3528266,
half dt8.3456086. Half-dt norm ratio.99915434 (-.08457%); half-amplitude
ratio1.00001850. Corresponding field L2 differences are1.73529% and1.45573%.
Growth remains under both refinements. This does not prove convergence of
the full field, a formal eigenmode, or late-time physical accuracy.
Evidence: canonical_coupled_t2p8_convergence_2026-10-18.json.

## t3.0 refinement and completed t3.4 half-amplitude comparison

At t3.0004059, coarse departure squared growth16.4596074, half amplitude
16.4535463, half dt16.4061588. Half-dt squared norm ratio.99675274
(-.32473%), field L2 difference2.16663%; half-amplitude field difference
1.79241%. Scalar growth persists; spatial discrepancies remain.

Both half-amplitude continuations completed320 accepted finite/no-retry
steps and exited0 at t3.4004059. Squared growth79.6917085 versus coarse
79.6248769 (+.08393%); field L2 difference3.19823%. This rules against a
large-seed-only explanation of this diagnostic growth, not nonlinear or
non-smooth effects in general. Fine-dt continuation is still running.
Evidence: canonical_coupled_t3p0_convergence_2026-10-18.json,
canonical_coupled_t3p4_amplitude_2026-10-18.json.

## t3.2004059 intermediate comparison completed

Matching amplitude checkpoint663 and fine checkpoint1142 were analyzed.
Departure squared growth is30.66059768 for coarse,30.55895536 for half
amplitude, and30.59984406 for half dt. Relative to coarse, squared norm
ratios are0.99668492 and0.99801851, respectively. Departure field relative
L2 differences remain1.93771% and2.45428%. Thus the roughly30-fold squared
response growth persists under these refinements; a timestep artifact alone
does not explain it. This is not a formal convergence order or unstable
eigenmode proof. Different restart histories remain a comparison caveat.
The fine-dt final t3.4 endpoint remains pending; historical t6-10 and long
physical validation remain required.
Evidence: canonical_coupled_t3p2_convergence_2026-10-18.json.

## t3.4004059 full refinement completed

All five branches exited0. Final analyzer verifies320 accepted steps per
amplitude branch and640 per fine-dt branch, all finite with zero retries,
matching physical times and persistent IDs/masses. Squared departure growth
is79.62487689(coarse),79.69170850(half amplitude),80.19384733(half dt).
The fine/coarse squared norm ratio is1.00714564 (+0.71456%), but departure
field relative L2 difference increases to7.32003%, from2.45428% at t3.2.
Half-amplitude departure field difference is3.19823%. Fine/coarse position,
velocity and pressure field differences are0.89395%,1.68422%,1.49512%.

The strong growth survives refinement. However, agreement of aggregate
norms does not establish spatial response convergence; the growing7.32%
discrepancy must remain explicit in mechanism claims. No formal convergence
order or unstable eigenmode is established, and restart histories differ.
The hard reconstruction cutoff already has independently demonstrated local
force discontinuities, but this comparison does not attribute the entire
field discrepancy to that cutoff. Retain the running smooth/legacy physical
comparison through t10 as the next decision gate; do not extend these short
response diagnostics merely to replace historical physical failure tests.
Evidence: canonical_coupled_t3p4_convergence_2026-10-18.json.
