# MPI design for bounded Kepler wall events

Date: 2026-10-31

## Scope

This is a code-specific design audit for making the existing bounded-W0 RK4
wall-event transaction compatible with more than one MPI rank. It does not
change the active executable or claim that the MPI path is implemented. The
current long-time W0 run remains the single-rank reference until an isolated
candidate passes the gates below.

## Current source contract

`Exam/gfs_wall_events.inc:gfs_wall_begin` rejects `NID(s) != 1`. The reason is
broader than the event-time loop: the transaction snapshots owned particle
payloads, builds an ID-indexed accepted-state table, scans RK4 dense paths for
the earliest polygon crossing, retries at a refined event time, applies slip
impulses, and closes global momentum, angular-momentum, and energy ledgers.
The current helpers treat all owned particles as one globally complete set.

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

The particle partition must remain fixed throughout one RK4 trial, or the
snapshot/identity transaction must be redesigned to follow stable IDs across
rank migration. Verify this from the actual sorting and migration call sites;
do not assume it from the fixed local `VORO_NP` check.

## Candidate implementation order

1. Keep the one-rank path byte-for-byte behaviorally unchanged by placing MPI
   collectives behind a separate, opt-in `GFS_KEPLER_WALL_EVENTS_MPI=1`
   contract. Reject inconsistent environment settings across ranks.
2. Factor local scans from global decisions. Keep snapshots and geometry
   caches rank-local; reduce only compact event and conservation records.
3. Add deterministic minimum-event selection and broadcasted root state.
   Fail closed on duplicate IDs, ties that cannot be resolved, rank migration,
   or a missing event owner.
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
