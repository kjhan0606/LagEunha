# Opt-in OpenMP validation for bounded-wall transactions

The active Q1 thermal64 wall-event candidate checks every particle's ID,
fixed mass, W0, finite state and thermal payload repeatedly across RK4
stages. Its baseline `gfs_wall_valid` performs those checks serially. A
separate [candidate patch](KEPLER_WALL_GUARD_OMP_2026-10-29.patch) adds a
default-off `GFS_WALL_OMP_VALIDATE=1` path. It retains the serial path when
the option is absent or OpenMP is unavailable. The parallel path claims each
ID with an atomic generation exchange, rejects duplicate/missing/out-of-range
IDs and changed mass, and reduces invalid-state flags before any trial can
be accepted. If multiple different failures occur concurrently, the chosen
diagnostic reason can differ from serial ordering, while both reject.

The patch is against the isolated thermal64 candidate's
`Exam/gfs_wall_events.inc` with SHA256
`7a831bfac2a6ac4c1823f253309f1f0d31835899ea76fbb2e95c139c03138be4`.
Candidate header hash is
`05b1bd283812c9de867c4ae99bf99fb925dd2f9a5ae4a92e8f29ebb9a08ebd88`.
The CUDA+OpenMP binary has SHA256
`da70844984a52fec77e2d6dec451e9b5695e761ba93a42d0f0d079c1091edf44`;
the preserved GPU object matches the original thermal64 candidate exactly.
An initial erroneous build lacked `USE_CUDA` and was replaced before any run.

Using the same binary and saved Q1 thermal64 starting checkpoint on Ada,
the no-impact four-step serial/OMP comparison passes byte-for-byte complete
native records except runtime next pointers, complete footer and printed RK4
and wall ledgers. A separate one-step forced-rejection pair passes the same
comparison and verifies exact rollback and device-cache invalidation in both
branches. See [four-step receipt](KEPLER_WALL_GUARD_OMP_PAIR_2026-10-29.json),
[retry receipt](KEPLER_WALL_GUARD_OMP_RETRY_2026-10-29.json), and
[auditor](audit_wall_guard_omp_pair_20261029.py). The executable contains the
OpenMP option and each corresponding frozen runner sets it explicitly.

The [extracted-function harness](test_wall_guard_omp_20261029.py) ran both
serial and parallel paths with 2 and 8 OpenMP threads; it rejects 12 separate
synthetic faults including duplicate/out-of-range IDs, changed mass, invalid
thermal state and invalid geometry. This is a mock harness using the actual
extracted guard function, not native fault injection or race-detection proof.

Both native comparisons ran concurrently with the ongoing Ada production
diagnostic, so their timings are **not valid speed measurements**. Multiple
realistic checkpoint geometries, an exact
restart bridge, and isolated ABBA performance test remain required before
the option is promoted. The current wall-event branch itself is one MPI rank;
this change does not establish MPI scaling, spin closure or Kepler stability.
The running Q1 binary was not changed or restarted.

## Actual thermal64 wall impact

Using the separately preserved compact-ID Q2 checkpoint at t=48.37534965270793,
the same candidate was run with serial and OpenMP guards and thermal64 enabled.
The native initialization promotes the stored float internal energy exactly
and recalculates double volume before stepping; no separate initialization
checkpoint was written. Both runs reach t=48.378, hit the same facet, and
report eight rollbacks and four accepted prefixes. The physical, finite-state,
W0/fixed-mass, double-U/V payload and EOS gates pass. Total-energy change
from the stored pre-initialization state using the authoritative final U64 is
-4.61e-14 relative in both cases. All records except runtime next pointers,
and the full footer, agree between serial and OpenMP. See
[thermal collision receipt](KEPLER_WALL_GUARD_OMP_THERMAL_IMPACT_2026-10-29.json)
and its [auditor](audit_wall_guard_omp_thermal_impact_20261029.py).

The serial candidate endpoint also matches the pre-existing thermal64 wall
candidate's independently saved t=48.378 endpoint under the same full-record
comparison; see [original-candidate parity](KEPLER_WALL_GUARD_OMP_THERMAL_ORIGINAL_PARITY_2026-10-29.json).
This adds a physical-impact check in thermal64, but the geometry is one
preserved Q2 event, not a Q1 boundary-convergence or long-time test. These
Ada runs shared the device with the ongoing Q1 simulation and cannot qualify
performance. Exact restart-bridge, independent ABBA timing and MPI remain open.

## Actual wall-impact comparison

The saved Q2 gapped-ID t=48.37534965270793 checkpoint was replayed to t=48.4
with the same candidate binary, serial and OpenMP validation. Both hit the
same wall face at t=48.37702401924956, executed 13 real rollbacks and 17
accepted prefixes, passed the existing physical gates, and produced identical
full restart records excluding only runtime next pointers and identical full
precision footers. See [impact receipt](KEPLER_WALL_GUARD_OMP_IMPACT_2026-10-29.json).
This collision uses the legacy thermal representation of that checkpoint;
the earlier Q1 no-impact and forced-retry checks exercise thermal64. On its
own it does not establish a thermal64 physical-impact replay, MPI, a restart bridge, or
long-time accuracy. The GPU was shared and timing is not a speed qualification.

In the raw serial-case log, an unrelated CPU line was appended immediately
after an `sfl_cum` numeric field. The original auditor stopped when parsing
that field. A separate [interleaving-aware auditor](audit_kepler_wall_events_interleaved_20261029.py)
accepts a complete numeric token followed by a new log marker and rejects a
truncated number; the unchanged raw logs were then audited using
[this wrapper](audit_kepler_wall_guard_omp_impact_20261029.py). Tests cover
the interleaved zero, a nonzero value and a truncated exponent. The original
auditor was restored byte-for-byte because pending H200 job 440010 freezes
its hash; that job does not use this modified copy. The old positive-control
result remains a separate existing result.
