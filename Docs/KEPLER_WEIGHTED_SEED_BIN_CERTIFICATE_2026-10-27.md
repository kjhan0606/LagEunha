# Conservative seed-bin exclusion for W0 and weighted cells

The previous [work-count diagnostic](KEPLER_Q2_GPU_TESS_WORK_DIAG_V8_2026-10-27.md)
suggested skipping bins whose particles cannot cut the current polygon.
The simple `d_min > 2R` certificate applies to W0; a power diagram needs
the weight difference as well. The following is a proposed exact exclusion,
not a promoted production optimization.

For an owner of weight `wi`, a candidate of weight `wj`, separation `d`,
and the existing clamped fractional plane position,

```
f(d,wi,wj) = clamp(0.5 + 0.5*(wi-wj)/d^2, 0.05, 0.95)
plane_distance = f*d.
```

Let `Delta=max(0, global_max_weight-wi)`. If the global bound includes
every real and padded candidate, `wi-wj >= -Delta`. For `d>0` the lower
bound

```
L(d) = max(0.05*d, min(0.95*d, 0.5*(d-Delta/d)))
```

is nondecreasing. Let `d_min` be a conservative lower bound on the distance
from the owner to **every actual point** in a bin, and let `R` be the
current maximum polygon vertex radius. When `L(d_min)>R`, every member
already fails the existing radius rejection, so the bin can be omitted
without changing the sequence of actual clipping operations. The later
adaptive-completeness search must remain unchanged. For zero weights this
reduces to `d_min>2R`.

The standalone host/CUDA-compatible
[certificate](gfs_seed_bin_certificate_20261027.h) uses squared distances,
conservative coordinate uncertainty, strict comparison, and treats edge
bins as unbounded in their edge coordinate. A small independent randomized
[host test](test_gfs_seed_bin_certificate_20261027.cpp) passed 500,000
weighted and edge-bin cases, certifying 65,040 skips with no observed unsafe
skip. This test samples only a tiny subset of inputs and does not prove the
native data structure satisfies the certificate's assumptions.

The native linked-list code computes a linear index from separate integer
coordinates. A linear index can in principle fall within the allocated
array even if one coordinate is out of range, so it is **not sufficient**
to assume that every CSR entry belongs to the bin's geometric rectangle.
Before integrating the exclusion, validate the actual CSR positions against
their assigned bins on every stage, including padded particles and MPI halo
states. If any entry or weight bound fails, the implementation must fall
back to the complete existing enumeration; that validation must be cheap
and race-free. Check floating-point margins and the global weight maximum
for all ranks, and retain the existing geometry aborts.

Only after those invariants are enforced should a candidate be built and
tested against the existing implementation on identical W0 and weighted
checkpoints, including distorted/jittered states and boundary cells. Compare
all fields, geometry and actual stage time on the same GPU; do not infer a
speedup from fewer candidate visits. This optimization has not been applied
to any running simulation.

The isolated [CUDA patch](apply_kepler_weighted_seed_certificate_20261027.py)
now adds an opt-in `GFS_GPU_SEED_BIN_CERT` candidate to the v5 tree. A
stream-ordered GPU kernel first checks every CSR entry's index, geometric
bin membership with a conservative margin, and weight against the global
upper bound. Any invalid entry sets a device flag that causes all owner
threads to retain the original complete seed scan for that stage. The
adaptive search remains unchanged. CPU Slurm **437511** builds this separate
candidate and reruns the 500,000-case host test; build success is not proof
of device fallback, field equivalence, or useful GPU speedup. Running jobs
retain their existing binaries.

The instrumented follow-on build **437513** completed `0:0`. It adds
opt-in per-stage telemetry for invalid-flag state and the number of bins
skipped, plus a diagnostic forced-fallback flag. A10 job **437517** was
submitted using the [pair stager](stage_kepler_seed_certificate_pair_20261027.py).
It compares baseline, enabled-certificate and forced-fallback replays from
the same jittered W0 `t=6.5` checkpoint and the same weighted `t=0.825`
checkpoint. Each native replay has a physics audit; the
[pair auditor](audit_kepler_seed_certificate_pair_20261027.py) requires
ten fields to match bitwise, records actual device skip counts, and checks
that forced fallback disables skips. Forced fallback is a flag-route test,
not an injection of a genuinely corrupt CSR. Diagnostic atomics affect
timing, so production speed must be measured separately if correctness holds.

GPU job **437517** completed `0:0`. All six native replay physics audits
passed, and the [pair receipt](KEPLER_SEED_BIN_CERTIFICATE_GPU_2026-10-27/pair_audit.json)
confirms bitwise equality of the ten compared fields for both enabled and
forced-fallback paths in W0 and weighted cases. The actual device telemetry
reports approximately `1.195 million` skipped bin visits per RK4 stage in
each enabled case and zero skipped visits in each forced-fallback stage.
This does not establish that those bins were populated, that an actually
malformed CSR is handled correctly, or that the candidate is faster. The
instrumented run times were only one step per variant and include diagnostic
atomics and startup. A separate same-GPU uninstrumented timing test remains
necessary, as do later/distorted states and MPI halo validation before any
production promotion.

The separate uninstrumented same-A10 job **437520** completed `0:0`.
The [speed audit](KEPLER_SEED_BIN_CERTIFICATE_GPU_2026-10-27/speed_audit.json)
passes all native gates and compares ten fields bitwise equal after 13 RK4
steps per variant. Excluding the first four of 52 tessellation calls, the
median host-and-device tessellation time was:

| Case | Existing scan | Certificate | Relative change |
|---|---:|---:|---:|
| W0 jitter | 162.090 ms | 165.326 ms | +2.0% |
| Weighted jitter | 172.146 ms | 175.178 ms | +1.8% |

Whole-run elapsed times were `23.982→24.440 s` and `24.850→25.602 s`,
respectively. These are short sequential measurements with the baseline
first, not broad hardware or long-time benchmarks. Despite demonstrating
actual bin skips and field equality, the candidate did not improve this
measurement and is **not selected for the running calculations**. Fewer bin
visits alone did not translate into a speedup; the current experiment does
not isolate whether empty bins, validation cost, extra arithmetic or another
factor explains the difference.
