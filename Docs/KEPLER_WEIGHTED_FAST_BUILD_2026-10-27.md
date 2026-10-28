# Weighted capacity64: isolated acceleration candidates

The validated bounded W0 implementation includes host-side timing and
configuration-cache improvements, with bitwise paired tests described in
[the fast RK4 record](KEPLER_BOUNDED_FAST_RK4_2026-10-27/README.md).
Those W0 tests alone do not validate the weighted code.

CPU build job **437501** tried all four W0 timing/acceleration patches on
the isolated v5 weighted source with `patch --fuzz=0`. The first timing
patch failed four of 17 hunks: the weighted stage code has diverged. The
job stopped before compilation; no running binary or source checkout was
changed. Its `patch_performance.log` in the isolated v6 build directory
preserves the failed dry run. It is not safe to force these patches.

CPU build job **437503** tries the independent Kepler acceleration parameter
cache patch alone on the identical v5 weighted inputs. Its staging runner is
[stage_kepler_weighted_accel_cache_v7_20261027.py](stage_kepler_weighted_accel_cache_v7_20261027.py).
Even if compilation succeeds, the candidate must be compared against v5
from the same accepted weighted checkpoint, including native field equality,
physics gates and matched runtime, before any longer weighted calculation
uses it. No speed or physical-equivalence claim for the weighted path has
been made yet.

Build **437503** completed `0:0`; the cache patch applied without fuzz,
the projection unit tests passed, and its binary SHA-256 is
`52da3334cd9452eaef3371253c10c7cebaf2ca2683e7e2e5a41ac3f5e256ae54`.
The cache initializes process-static data under the current serial call
assumption; a future OpenMP parallelization of this call would require
explicitly thread-safe initialization.

The [pair stager](stage_kepler_weighted_cache_pair_20261027.py) prepared
A10 job **437506** with `afterok:437497`. It runs the existing v5 and the
candidate consecutively on the same verified `t=0.0456165` checkpoint to
exact `t=0.06`. Each run has its own 21-check native physics audit; the
[pair audit](audit_kepler_weighted_fast_pair_20261027.py) additionally
compares ten physical fields bitwise at the initial, parent and final states.
The calculation is a short controlled equivalence check; its run times are
not sufficient to claim production performance or long-time equivalence.

Job **437506** completed `0:0`. Both 21-check physics audits passed and
the [pair audit](KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/cache_pair_t006.json)
found all ten compared fields bitwise equal at the initial, parent and final
states. The whole-run wall times in this short sequential pair were 16.642 s
for v5 and 15.320 s for the cache candidate. They include startup and I/O,
so the 8% difference is not a sustained performance result. The running
v5 chain is unchanged; the experimental candidate has not been promoted to
a long simulation.

As an alternative to the [unsuccessful seed-bin optimization](KEPLER_WEIGHTED_SEED_BIN_CERTIFICATE_2026-10-27.md),
the existing W0-tested wall-select GPU patch was applied with zero fuzz on
top of the cache candidate. CPU build **437525** completed `0:0` for this
independent v10 candidate. Its [stager](stage_kepler_weighted_wallselect_cache_v10_20261027.py)
and the two independent A10 comparisons are:

- **437528** compares v5 and v10 from the weighted accepted step-start at
  `t=0.0456165` through `t=0.06`, crossing the earlier capacity-failure area.
- **437527** measures v5 and v10 consecutively on the same GPU over the
  jittered W0 `t=6.5→6.525` and weighted `t=0.825→0.85` intervals; both
  physical gates and all ten compared fields are required to pass.

The original source patch has prior W0 evidence, but these weighted tests
are needed because weighted cell geometry and force work can respond
differently to a change in clipping order. Running v5 jobs are unmodified.

Job **437527** completed `0:0`; its
[matched-speed receipt](KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/wallselect_cache_speed_audit.json)
passes all native gates and finds all ten compared fields bitwise equal at
the end of the 13-step replays. Excluding the first four of 52 tessellation
calls, median time was `163.007→149.731 ms` for W0 (−8.1%) and
`172.687→159.772 ms` for weighted (−7.5%). The complete sequential run
times were `23.852→20.430 s` and `25.185→21.579 s` respectively. These
measurements compare v5 with the combined cache and wall-select changes;
they do not independently attribute each contribution. The baseline ran
first, and the measurements cover short one-GPU checkpoints only. Job
**437528** separately checks the earlier weighted capacity-failure point;
neither pair establishes late-time weighted or MPI-rank equivalence.

Job **437528** also completed `0:0`; its
[capacity-point pair audit](KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/wallselect_cache_capacity_audit.json)
passes both physics gates and finds all ten compared fields bitwise equal
at the weighted `t=0.0456165→0.06` replay. Its whole-run times were
`16.590→14.257 s`, again a short sequential observation.

The [late-pair stager](stage_kepler_weighted_wallselect_late_pair_20261027.py)
submitted **437538** with `afterok:437499:437527:437528`. It verified the
accepted weighted `t=6.5` checkpoint, earlier pair receipts, both binary
hashes and input lineage, then compared v5 and v10 on the same A10 through
`t=6.525`. Job 437538 completed `0:0`; both native audits passed **21/21**,
and the pair receipt found all ten fields bitwise equal at initial, parent
and endpoint states. The [late-pair receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_wallselect_weighted_late_pair_v1/pair_audit.json)
records the comparison. Whole-run wall time was `24.159→20.867 s` (13.6%
lower for v10) on this short sequential replay. This supports using v10 as a
performance candidate; it does not establish long-time equivalence or MPI
rank equivalence. The accepted v5 `t=10→20` run remains unchanged. Its
pending H200 `t=20→30` continuation is also left intact to measure the
baseline binary's long-interval H200 throughput before a separately gated
v10 continuation.

An isolated **v11** build layers the diagnostic-only thermal64 flag hoist
onto v10. The patch applies to the exact v10 source archive with zero fuzz;
the geometry projection unit test passes. The binary SHA-256 is
`adeca5d738bce85a740b966fbbef04e3aaaba898ac7978eb5db797a620eaf5cd` and
the modified `Exam/exam.c` SHA-256 is
`9c66f067401bb866c9fe04b3d39b02fda5029b921c24a5d5eefb0454d1d534f7`.
Build provenance is in the [v11 receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_wallselect_thermalhoist_capacity64_build_v11/build_receipt.json).
It has not yet passed a native weighted pair. H200 job **437623** will
compare v10 and v11 from the accepted v5 `t=30` state to `t=30.025`; its
dependent **437626** continues v11 to `t=40` only if the pair gate succeeds.
