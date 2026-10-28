# Conditional A10 sparse Kepler continuation beyond t=10

The corrected narrow and wide CUDA+MPI+OpenMP W=0, fixed-mass runs are
currently evolving to `t=10` on separate A10 cards (jobs **437254** and
**437253**). Their accepted binary includes `sm_86`; no CPU or Python gas
solver is substituted. The [first native intermediate comparison](KEPLER_SPARSE_A10_T165024_INTERMEDIATE_2026-10-27.md)
passes disk state and early orbital checks, but a separate radial audit
already measures about `0.2653` absolute density change near the outer
periodic box faces. Whether that disturbance reaches the disk remains
unresolved.

## Measured capacity, not a run-time promise

At the 2026-10-27 05:48:57 KST Slurm check, wide reached step 984,
`t=1.90454`, after 16:40 elapsed; narrow reached step 961,
`t=1.85989`, after 14:28. Both retained `dt≈0.00193984`. In the recent
91-second interval both advanced 41 steps, about **2.22 wall seconds per
accepted step** with both jobs sharing `syn02`. Holding that speed and
step size would put `10→60` near 15.9 hours and `60→120` near 19.1 hours.
The step size, contention, halo cost, and physical state may change, so
these are provisional capacity estimates. A single `10→120` interval has
less margin against the 48-hour Slurm cap; a `t=60` checkpoint is the
safer first long gate if `t=10` passes.

The `a10` partition's `MaxTime=2-00:00:00` was checked directly with
`scontrol show partition a10 -o`; the runner's 48-hour request is valid.
The earlier sparse H200 `.001→.825` path averaged about `0.465 s/step`
with halo 12, while the concurrent A10 jobs measured about `2.22 s/step`
with halo 32 near `t≈1.9`. Different halo widths, epochs, and node
contention prevent using their ratio as a hardware speedup factor. Choose
the long-run device after the actual `t=10` throughput and H200 quota are
known; do not duplicate the same long evolution on both devices.

The [offline local-halo reconstruction](KEPLER_Q2_SPARSE_LOCAL_HALO_FALLBACK_T2221_RECONSTRUCTION_2026-10-27.json)
at exact-orbit `t=2.221` used 9,091 periodic images versus 53,074 for a
global halo 32 and reproduced the affected cell areas to relative
`1.98e-13`. This is a promising performance target but only a geometry
experiment. It does not validate the native GPU face, force, work or halo
state path, so the accepted long runner retains halo 32.

## Prepared A10 bundle

The isolated [A10 runner](kepler_sparse_area_long_a10_20261027.slurm)
and [stager](stage_kepler_sparse_long_a10_20261027.py) derive from the
existing H200 long-run contract. A source diff verifies that only the
job name, A10 partition/GRES, device log, and read-only exact-orbit
Voronoi and seam-provenance diagnostics were added. EOS, force, RK4, timestep, halo 32,
fixed mass, W=0, initial state, and the native executable stay the same.
Both prospective long runners now set `GFS_GPU_TESS_PROFILE=1` to print
host/device tessellation timing for the first eight calls and every
1,024th call; this is diagnostic timing, not a physics switch.
The runner refuses a predecessor unless both its physical and independent
state receipts pass, then verifies the checkpoint, binary, and parameter
hashes before copying a restart. Each endpoint must pass the native
orbit/energy/entropy gate, the state consistency gate, the
[native-versus-orbit density diagnostic](audit_kepler_sparse_native_ballistic_20261027.py),
and the [outer-source provenance diagnostic](audit_kepler_seam_provenance_20261027.py).

An isolated scripts-only bundle at
`/tmp/LagForce_sparse_long_a10_stage_20261027_0554` passed Python import,
`bash -n`, and every staged `sha256sum -c` check before the profiling flag
was added. Stage a new immutable bundle after choosing the device; the
current A10 runner SHA-256 is
`1ccc6ecafac8be71ae016ea48c612f7284ef10daa6bc11525c06dc987fa42db7`.
**No continuation beyond `t=10` has been submitted.** The predecessor's
complete `t=10` checkpoint and receipts do not yet exist.

When `t=10` completes, first check its native physics receipt, run the
independent state and exact-orbit density audits on the final snapshot,
and inspect the radial seam-to-disk trend. CPU jobs **437262–437265** now
queue the `t≈6.6` and `t=10` radial/geometry comparisons from the existing
A10 trajectories. If these pass and the boundary trend is acceptable,
stage a fresh
`10→60` A10 work directory and submit with a realistic wall cap and
intermediate dumps. Promote `60→120`, then later segmented `t=600` and
the explicitly defined longest-orbit endpoint only after each predecessor
gate. A strong boundary-driven disk error requires a boundary-control
experiment before interpreting a long run as the GIZMO Kepler test.

## Gate halted by the GPU corner-cap failure

The two original A10 jobs **437253/437254** have now failed near `t=2.5`
at `[GPU_CLIP_CAPACITY] nc=48 capacity=48`; the independent last-checkpoint
radial and Voronoi area audits pass. See the [failure and checkpoint
record](KEPLER_SPARSE_A10_T10_PARALLEL_2026-10-27.md#both-a10-runs-stopped-at-a-gpu-clipping-capacity-2026-10-27).
No `t=10` predecessor exists, so do not stage or submit `10→60` yet.
Repair and validate the CUDA clipping path, then resume from the preserved
accepted checkpoints around `t=2.475`, without replaying the original
`.825→2.475` interval.
