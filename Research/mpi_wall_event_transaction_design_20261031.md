# MPI design for bounded Kepler wall events

Date: 2026-10-31

## Scope

This is a code-specific design audit for making the existing bounded-W0 RK4
wall-event transaction compatible with more than one MPI rank. It does not
change the active executable or claim that the MPI path is implemented. The
current long-time W0 run remains the single-rank reference until an isolated
candidate passes the gates below.

## Current source contract

`Exam/gfs_wall_events.inc:gfs_wall_begin` permits multiple ranks only with
`GFS_KEPLER_WALL_EVENTS_MPI=1`. The ordinary path rejects a multi-rank run
without that explicit opt-in. The event transaction snapshots each rank's
owned particle payload, scans RK4 dense paths, retries at a refined event time,
applies slip impulses, and closes global momentum, angular-momentum, and
energy ledgers.

The existing `GFS_WALL_PARALLEL_VALIDATION` option can parallelize one local
state scan, but it does not make the event transaction MPI-safe. The
diagnostic-only `GFS_WALL_RHS_PROBE_ONLY` also has an independent one-rank
guard. Do not remove either restriction as a shortcut.

## Required rank-wide decisions

Every rank must make the same accept/reject and retry decision at every RK4
trial. Add collectives at these boundaries:

1. **Stage validity:** reduce local geometry, band, identity, and finite-state
   failures with a communicator-wide OR before any rank enters the next stage.
2. **First wall event:** reduce local dense-path crossing fractions to one
   global minimum, then select a deterministic `(owner rank, stable ID, side)`
   for ties. Broadcast that tuple so all ranks enter the same root-refinement
   state.
3. **Root refinement:** the owner computes the selected particle's signed
   gap and normal speed; broadcast both before updating the common bracket.
   Recompute global competing crossings on each trial. A rank-local root or
   competing event must never cause a rank-local return.
4. **Contact validation and impulse:** collectively reject unsupported
   outside endpoints or corners. Apply each resolved contact impulse on the
   rank that owns that particle exactly once; reduce collision impulse,
   torque, roundoff, and impact count before acceptance.
5. **Conservation ledgers:** `gfs_wall_totals`, pressure/gravity impulses,
   and commit residuals must use global sums. Only rank zero should emit the
   aggregate `[WALL_ACCEPT]` line.
6. **Rollback:** restore each rank's local accepted payload and local saved
   ledgers, then collectively verify that all ranks restored the same driver
   time, timestep, global owned count, and transaction decision before retry.

The particle partition changes during RK4 predictor calls to
`migrateTreeVorork4Particles()`. Rank-local snapshots remain valid for rollback
because `gfs_ssp_restore` reinstates each rank's saved pre-trial payload and
count. The MPI candidate now allows accepted ownership changes, builds a
replicated immutable global-ID/mass map once, checks local particle validity
and global cell count at each stage, and reconstructs accepted start positions
from RK4 weighted displacements when scanning the final dense path.
`GFS_WALL_MPI_VERIFY_IDS=1` adds an exact global ID-set reduction immediately
after each fixed-polygon migration; it is intended for short qualification
runs because it communicates a dense ID vector per migration. These are
source changes only; ownership migration has not yet passed a runtime test.

## Candidate implementation order

1. Keep the one-rank path byte-for-byte behaviorally unchanged by placing MPI
   collectives behind a separate, opt-in `GFS_KEPLER_WALL_EVENTS_MPI=1`
   contract. Reject inconsistent environment settings across ranks.
2. Factor local scans from global decisions. Keep snapshots and geometry
   caches rank-local; reduce only compact event and conservation records.
3. Add deterministic minimum-event selection and broadcasted root state.
   Keep failing closed on wall contacts until this path and the distributed
   impulse/ledger behavior are qualified.
4. Add rank-wide fault injection for one retry. A single-rank forced rejection
   is insufficient: every rank must restore and retry the identical interval.
5. Build a source-pinned binary and run from the same t=17.8 restart used by
   the fixed-mass reference. Compare one rank against two and four ranks at the
   same endpoint, recording density, velocity, energy, angular momentum,
   wall ledgers, accepted timestep sequence, and wall time.

## Release gates

- One-rank candidate reproduces the existing control before testing MPI.
- Two- and four-rank no-impact runs pass all current physics gates and agree
  with the one-rank endpoint within the established numerical tolerances.
- A deliberately rejected RK4 trial restores each rank exactly and retries
  globally; no rank advances time or commits an impulse early.
- A controlled wall-contact case resolves one global event, applies impulses
  once, and closes global linear and angular momentum ledgers.
- End-to-end timing demonstrates useful speedup after communication and
  collective costs. Keep the one-rank executable as the fallback unless the
  candidate passes every gate.

The separate multi-rank face diagnostic is evidence about GPU face exchange,
not a substitute for these transaction and conservation gates.

## No-impact MPI prototype in progress

The source worktree now has an opt-in no-impact prototype using
`GFS_KEPLER_WALL_EVENTS_MPI=1`. It synchronizes event-mode settings once,
keeps rollback snapshots rank-local, allows `pmigrate` to change local cell
counts, checks global cell count and replicated ID/mass membership, reduces
stage validity and endpoint violations, and globally sums start/end
conservation totals and wall ledgers. `GFS_WALL_MPI_VERIFY_IDS=1` additionally
checks the complete global ID set after each migration. If any rank predicts a
polygon-wall crossing, the trial still fails closed; the prototype does not
select a global owner, refine a distributed event root, or apply impacts across
ranks.

This candidate is only suitable for a short no-impact parity/performance
comparison. It is not a production wall-event implementation and must not be
used for a long segment until an MPI event owner/root/impulse path is added and
the controlled-contact gates pass. On 2026-10-31, the updated CUDA-enabled host
translation units `Exam/exam.c` and `Exam/KH/kh.c` compiled successfully with
the pinned Intel MPI/CUDA modules and no warnings. This was not a full
executable link or runtime migration/parity test; those gates remain pending
and the active baseline executable is unaffected.
