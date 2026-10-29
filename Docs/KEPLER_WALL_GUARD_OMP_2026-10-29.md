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
diagnostic, so their timings are **not valid speed measurements**. A physical
wall-impact comparison, multiple realistic checkpoint geometries, exact
restart bridge, and isolated ABBA performance test remain required before
the option is promoted. The current wall-event branch itself is one MPI rank;
this change does not establish MPI scaling, spin closure or Kepler stability.
The running Q1 binary was not changed or restarted.
