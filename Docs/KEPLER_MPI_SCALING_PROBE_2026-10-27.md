# Weighted Kepler MPI scaling probe — 2026-10-27

## Purpose

Measure the same weighted Kepler evolution with one and two MPI ranks on two
H200 GPUs, starting from the same independently accepted (t=40) checkpoint.
This is a short throughput and rank-consistency probe; it does not establish
long-time stability or a general MPI speedup by itself.

## Physics and comparison contract

- Binary: v11, SHA-256
  `adeca5d738bce85a740b966fbbef04e3aaaba898ac7978eb5db797a620eaf5cd`.
- Initial state: same original (t=0) snapshot and accepted weighted (t=40)
  restart for both branches.
- Rank counts: one-rank control and two-rank candidate; each rank uses eight
  OpenMP threads. The Slurm allocation requests `--gres=gpu:H200:2` and the
  batch script checks that exactly two GPUs are visible.
- Endpoint: (t=40.25). Each branch must pass all 21 native physical checks.
  The comparison records per-field bitwise status, maximum and RMS differences,
  scale-normalized RMS differences, and wall-time speedup. Rank-dependent
  floating-point differences are measured rather than required to be zero.
- The (t=40) checkpoint gate must pass first. The MPI job is chained with
  `afterok:437626`.

## Current state

- A10 job **437531** completed the exact `t=10→20` continuation with exit
  code 0. Its authoritative receipt is
  `/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_fresh_to20_a10_v1/audit_t20.json`.
  It passes all 21 configured checks: exact fixed mass, zero initial weights
  with evolved positive weights, disk density error `1.29209e-4`, orbit
  position error `1.17038e-6`, relative total-energy error `1.60557e-8`, and
  zero recorded injection. This validates W=0-to-Laguerre evolution for this
  finite interval, not a rotation model or local moving-face conservation.
- The separate fixed-polygon W0 job **437441** failed at `t=48.3753`, step
  11580, when site 55283 crossed the static polygon boundary by `3.15984e-6`.
  Its last completed RK4 ledger had relative total-energy drift `-1.108e-8`,
  no nonpositive internal-energy cells, and zero floor injection. This was a
  boundary-guard termination, not evidence of a disk energy blow-up. Its
  dependent `t=60→100` job **437442** was canceled before allocation because
  `afterok:437441` could not be satisfied.
- The full-domain periodic W0 job **437076** remains RUNNING on H200. A later
  log read reached `t=1.95547`, global step 26248; a nearby complete ledger at
  `t=1.95290` reported `dEtot/|Etot0|=-2.523e-10`, zero nonpositive
  internal-energy cells, and zero floor injection. The last 100 accepted-time
  markers had median `dt=7.642e-5`; several one-step minima recovered on the
  following steps, so there is no observed sustained timestep collapse. This
  is the relevant full-domain control for the historical `t>6` question; the
  halo-28 run has unchanged annulus resolution, as documented in
  [the resolution audit](KEPLER_HALO28_DISK_RESOLUTION_AUDIT_2026-10-27.md).
- Job **437613** completed the weighted `t=20→30` continuation; its t=30
  receipt passes all 21 configured checks. Pair job 437623 failed in
  restart staging before either native simulation, as detailed below.
- At `t=30`, the accepted weighted state has `w2_min=0` and
  `w2_max=1.28714e-5`; the disk's mass-weighted density change is
  `5.84257e-4` (gate `1e-3`), entropy-log change `2.29979e-5`, orbit-position
  error `4.16203e-6`, disk angular-momentum drift `-5.09213e-10`, and total
  energy drift `2.40220e-8`. This verifies finite-time variable-weight,
  fixed-mass evolution through t=30 for the accepted v5 branch; it does not
  yet establish t=600 survival or cell-spin behavior.
- The old dependent jobs 437626 and 437636 were canceled before allocation.
  Their fresh replacements are **437794** (matched pair), **437796** (v11 to
  t=40 after the pair gate), and **437798** (MPI rank comparison after t=40).
  The pair is currently pending H200 GPU quota; the later jobs wait on its
  acceptance gates.
- Duplicate job 437625 and the first MPI probe 437635 were canceled before
  allocation. MPI job 437636 was later canceled with its failed 437623
  dependency; replacement job 437798 is staged from the new t=40 gate.

## Interpretation limits

The current weighted model begins with every `w2` equal to zero, then evolves
nonnegative weights with at least one positive weight. Its fixed-mass contract
requires exact per-cell mass equality at the endpoint. Thus the initial
tessellation is Voronoi and subsequently becomes Laguerre. Passing the upcoming
MPI pair would show short-interval consistency and measured throughput for
this binary and machine configuration only. The
long-time GIZMO (t=600) goal, local moving-face conservation, and cell-spin
model remain unverified.

## Runner repair and replacement chain (2026-10-27)

The original pair job **437623** failed before entering either native run.
Its parent-restart Python heredoc used `root` without defining it; the log
ends with `NameError: name 'root' is not defined`, and Slurm records zero
elapsed time. The accepted t=30 parent receipt itself passes all 21 checks.
The pair runner now defines `root = str(Path.cwd())` inside Python, and the
MPI runner received the same correction. The failed v1 stage and its log are
preserved unchanged.

Fresh v2 stages were hash-checked, then submitted in order: **437794** is the
matched v10/v11 pair at t=30, **437796** is the v11 continuation to t=40 after
437794, and **437798** is the one-rank/two-rank H200 comparison after 437796.
The completed t=30 receipt was checked before 437794 was submitted without
an `afterok` dependency because Slurm rejected the already-completed 437613
dependency as a submission error. The old v1 dependent jobs 437626 and
437636 were canceled before allocation; their v2 replacements retain the
same physics and acceptance gates. No simulation replay or solver change was
made to repair the staging error.

## Matched-pair import recovery and t=40 continuation (2026-10-27)

The v10/v11 matched-pair job **437794** completed both native `t=30.025`
branches and each branch passed its 21/21 physics audit. Its final comparison
step failed only because the staged directory omitted the helper module
`audit_kepler_sparse_area_t0825_20261027.py` from Python's import path
(`ModuleNotFoundError`). The snapshots and receipts were preserved. Re-running
only the comparison audit with the `Docs` directory on `PYTHONPATH` produced
`pair_audit.json`; all four checks pass: both physics gates and bitwise equality
of the selected initial, parent, and endpoint fields. Endpoint `w2`, density,
and internal-energy maximum differences are exactly zero. This establishes
short-interval v10/v11 equivalence for those audited fields, not physical
closure or long-time accuracy. The paired runtime logs show 15.019 s for v10
and 11.451 s for v11 on this single `0.025`-unit interval (23.8% lower wall
time for v11). This short run includes startup and is only a preliminary speed
signal; the longer `t=30.025→40` continuation is needed to measure sustained
throughput.

The permanently blocked old jobs **437796** and **437798** were canceled
without allocation. The already-staged v11 continuation from `t=30.025` to
`t=40` was resubmitted as **437812**; its runner rechecks the recovered pair
receipt and candidate physics receipt before copying the restart. The matched
one-rank/two-rank probe was resubmitted as **437813** with `afterok:437812`.
At 21:05 KST H200 resources became available and replacement job **437812**
started. At the 21:11:32 log sample it had advanced from the accepted
`t=30.025` restart to `t=31.1714` (step 597, `dt=0.00193939`). The RK4 ledger
reports relative total-energy drift `-9.007e-10`, `dEhyd/Ehyd0=-2.097e-6`, and
zero floor/source-floor injection. The `t=40` endpoint audit is still pending.
Across 2,454 logged GPU weight-solve stages, the median was 4.881 ms, p95
6.493 ms, p99 8.039 ms; 15 stages exceeded 0.1 s, with a maximum of 0.381 s.
The MPI comparison **437813** remains gated behind 437812. No solver run was
repeated to repair the missing import.
