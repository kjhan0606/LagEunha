# Bitwise-validated bounded Kepler RK4 acceleration

This directory preserves the final `Exam/exam.c` overlay for the fixed
1,024-gon W=0 native CUDA+MPI+OpenMP Kepler branch. The file SHA-256 is
`6b1d54c9de7ed7ce710293d92daac7e055434acf09a3674dd4026b9cde0c58f9`.
All other source files remain those of the
[near32 overlay](../KEPLER_FIXEDPOLY_NEAR32_SOURCE_2026-10-27/README.md).
The matching GPU object is byte-identical to the near32 reference:
`e0d30024a9afc5abab3965ba39d01d04ea5171b01d527b7f9b905c5afd2ccd6c`.
The built executable SHA-256 is
`d92b6cd3b8491678cdfafa7610608378db89ca9bd16585a52a51aaf26f3b3367`,
preserved in the immutable
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/rk_diag_hoist_pair_t085_a10_v1/profile/eunha2.exe`.

Starting with near32 `Exam/exam.c`, apply these patches in order:

1. [Exclusive RK4 gap timers](../KEPLER_RK4_GAP_PROFILE_2026-10-27.patch).
2. [One-time Kepler acceleration configuration](../KEPLER_ACCEL_CONFIG_CACHE_2026-10-27.patch).
3. [Ten-call diagnostic timer](../KEPLER_RK4_STAGE_DIAG_TIMER_2026-10-27.patch).
4. [Per-diagnostic thermal64 flag hoist](../KEPLER_RK4_DIAG_THERMAL_HOIST_2026-10-27.patch).

Applying all four patches to the preserved near32 source in a fresh local
directory reproduced the final source SHA-256 exactly. The build used
`USE_CUDA=1 GFS_VARPRESS_W0=1`, two login-node build threads, Intel MPI,
Intel OpenMP and CUDA 13.0.2. The [build script](../build_kepler_fixedpoly_gap_profile_20261027.sh)
restores the accepted GPU object before linking. The `GFS_RK4_GAP_PROFILE`
environment variable controls detailed timing output; the additional
`MPI_Wtime` calls occur even when it is unset, but no physics expression or
input field has been changed.

## Native paired evidence

- A10 **437330** compared the original near32 executable and the timing-only
  executable over the same accepted `.825→.85` interval. All ten parsed
  physical fields matched bitwise and both native physical gates passed.
- A10 **437331** compared the timing-only and acceleration-config-cache
  executables over the same interval. All fields matched bitwise. In that
  sequential pair, post-warmup median RK4 wall fell `1.68204→1.41845 s`
  (15.7%).
- A10 **437335** measured the ten full-state stage diagnostics at a median
  `0.274181 s` per step; adding only the diagnostic timer retained bitwise
  physical equality.
- A10 **437336** compared the cached executable with the thermal64 flag
  hoist. Both native gates and all parsed physical-field comparisons passed.
  On the same A10 allocation, the post-warmup median RK4 wall fell
  `1.467731→1.1830635 s` (19.4%); the ten unchanged full-state diagnostic
  checks fell `0.2787475→0.0075715 s`.

The runs and [timing report](../KEPLER_RK4_GAP_PROFILE_2026-10-27.md) give
short, matched speed measurements. The active `t=.825→6.5` and dependent
`6.5→10` physics jobs retain the original executable; the faster binary
must enter only through a separately audited accepted restart. No `t>6`
stability or GIZMO-long-time claim follows from these performance checks.

## Guarded sustained continuation queued

The [fast restart stager](../prepare_kepler_bounded_fast_followon_20261027.py)
and [A10 runner](../kepler_q2_bounded_fast_to20_a10_20261027.slurm) are
preserved separately from the active physics executable. A temporary
staging smoke test using an already accepted `.825` native state passed
the parent audit, checkpoint SHA/time, original and faster binary SHA,
and every staged-input SHA check. No simulation ran in that smoke test.

Slurm **437343** was queued with `afterok:437328`: first the unchanged
bounded chain must pass `t=6.5` as **437317**, then its unchanged-binary
follow-on must pass exact `t=10` as **437328**. Only then does the 437343
runner verify the full `t=10` parent physics receipt and copy its accepted
checkpoint into a new directory. It selects the exact faster executable
SHA-256 `d92b6cd3b8491678cdfafa7610608378db89ca9bd16585a52a51aaf26f3b3367`,
evolves to exact `t=20` with 2.5-unit checkpoints, and runs the same
17-check bounded native physics gate. Its A10 request is one GPU, two CPU
threads, 32 GB and five hours.

## Accepted bounded aligned t=20 endpoint

Job **437343** completed `0:0` in 1h49m59s, continuing from the accepted
`t=10` checkpoint. Its [exact t=20 audit](audit_t20.json) passes all 17
checks. The disk mass-weighted mean absolute density change is `2.28372e-4`,
the native-minus-exact-orbit Voronoi residual is `3.46573e-5`, the absolute
entropy-log change is `2.08748e-5`, and the total-energy relative change is
`1.61036e-8`. The endpoint SHA-256 is
`cf2907c33718a04d1ac6f5db288115afd1d9dfe818286ffe514a099bcb4211db`.
The native run's recorded wall time was 6,594.8 seconds for `t=10→20`,
including its checkpoint output. This is a sustained runtime observation;
the short matched tests above remain the controlled speed comparison.

Its dependent job **437348** started the `t=20→30` continuation after the
successful gate. These are aligned, bounded 55,351-site W0 diagnostics;
they do not establish jitter robustness, full 384² periodic behavior, or
the paper's long-time density result.

The completed step-3958 checkpoint from running job **437348** was
independently audited at `t=27.501213097602502`. The
[live receipt](live_step3958_audit.json) passes 19/19 checks, including
fixed mass and W0, exact executable hash, no recorded geometry abort through
the endpoint, zero injection, orbit and energy checks. Disk mass-weighted
mean absolute density change is `3.02664e-4`, exact-orbit Voronoi residual
`1.34682e-4`, mean absolute entropy-log change `2.86013e-5`, and relative
total-energy change `2.19759e-8`. The residual is larger than the
`3.46573e-5` observed at t=20; both values are within the current broad
`1e-3` gate, but this trend warrants continued long-time monitoring.

The [auditor](../audit_kepler_bounded_native_evolution_20261027.py) now has
optional `--through-endpoint` log truncation and binary-hash verification;
its default final-endpoint behavior is unchanged. The live receipt uses the
log only through the saved step marker, and says nothing about later stages.
Re-running the revised auditor without the new options on the preserved
`t=20` state reproduced the existing receipt exactly, including all fields
and checks.

A sign-resolved decomposition from the same native snapshots gives mean
disk `Delta log K=-1.55115e-5` at [t=20](entropy_sign_t20.json) and
`-2.08347e-5` at [t=27.5012](live_step3958_entropy_sign.json). At the later
time, the positive and negative components are `3.88330e-6` and
`2.47180e-5` respectively; the inner `0.5≤r0<0.75` band has signed mean
`-1.01134e-4`. The signed whole-disk change is therefore negative in this
aligned control even though some cells have positive changes. These are
state differences, not an isolated irreversible-heating or face-work
ledger, and should not be used to diagnose the jittered disk.

Job **437348** subsequently completed `0:0` at exact `t=30`; its
[native receipt](audit_t30.json) passes all 17 configured checks. The
independent [endpoint-limited receipt](audit_t30_independent.json) also
passes the two additional geometry-abort and executable-hash checks. The
disk mean absolute density change is `3.66723e-4`, native-minus-ballistic
Voronoi residual is `2.10251e-4`, entropy-log absolute mean is `3.14251e-5`,
and relative total-energy change is `2.39930e-8`, with zero recorded energy
injection. The endpoint SHA-256 is
`12ea9ff96ec70c4816f834618229c972d0a8e47814bcdc2a38a20a8e89aad50e`.
This is the aligned, bounded W0 control, not the jittered weighted run or
the GIZMO `t=600` benchmark. The conditional `t=30→60` job **437441** is
has started from the audited t=30 checkpoint.
The [sign-resolved t=30 entropy state](entropy_sign_t30.json) has whole-disk
mean `Delta log K=-2.32438e-5`; 28.3% of disk mass has positive change.
In the inner `0.5≤r0<0.75` band the mean is `-1.06265e-4`. As above, these
are signed state changes, not a direct heating-rate or dissipation ledger.
