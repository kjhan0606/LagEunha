# Sparse buffered-area Q2 checkpoint at t=2.4754

## Why this branch matters

The 54,312-site sparse buffered-area run retains the same 49,038 cells in the
active `0.5<=r<2` disk as the full 384² Q2 reference, while reducing the
surrounding population. At exact `t=0.825`, its disk density change was nearly
identical to the full run, while its absolute entropy-log change was much
smaller. The `t=0.825` strict gate passed. This is a sparse-buffer control,
not the separate locally refined inner-resolution IC.

## Failure and restart point

The older H200 and A10 continuations both stopped at `t=2.50986`. Their logs
end with `[GPU_CLIP_CAPACITY] nc=48 capacity=48`,
`[GPU_FINAL_VOLUME] invalid=1`, then `EXIT:1`. The last complete snapshot is
`khout.001318.dat`, exactly `t=2.4753977144699952`; do not replay the earlier
trajectory.

The independent [checkpoint audit](KEPLER_Q2_SPARSE_PROFILE_T2P475_CHECKPOINT_AUDIT_2026-10-28.json)
compares it with the exact `t=0` state. The audited state has 54,312 complete
IDs, bitwise fixed mass, exact W=0, finite positive density/volume/internal
energy, and a total periodic volume of `36.0000000478`. Relative whole-box
energy drift is `1.84e-8` and the last available RK4 ledger has zero floor,
dual-energy, and entropy-floor injection.

| Initial-radius cohort | Mean absolute density change | Mean absolute entropy-log change | Mass-weighted orbit-position error | Relative Lz change |
|---|---:|---:|---:|---:|
| `0.5–0.75` | `4.996e-4` | `1.452e-5` | `7.34e-8` | `9.02e-13` |
| `0.75–1` | `2.399e-4` | `5.892e-6` | `3.54e-8` | `-3.77e-12` |
| `1–2` | `9.333e-5` | `3.806e-6` | `3.26e-8` | `2.17e-13` |
| `0.5–2` | `1.440e-4` | `4.922e-6` | `3.62e-8` | `-1.41e-13` |

The full-box Lz changes by `-0.9425%`, but a Cartesian periodic box is not
rotationally invariant. Use the disk cohort's angular momentum for this
diagnostic, not the full-box value as a conservation gate. The initial
`r=2.35–4.25` cells already show large edge motion (`1.26` mass-weighted
position error), so the good disk metrics at `t=2.4754` do not prove the
outer region is benign or that its later evolution cannot affect the disk.

## Continuation now queued

Job **438533** resumes from the audited checkpoint with the clip-swap candidate
binary (`6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0`),
uses H200:1, preserves the sparse-area W0/fixed-mass setup, crosses the old
`t=2.50986` failure, and targets exact `t=10`. It uses a fresh run directory
and runs the existing strict t=10 disk audit. At the 15:05:58 KST syntax
Slurm query, the job was pending for `QOSGrpGRES`; it had not started. Its
Slurm script is [here](kepler_q2_sparse_area_candidate_to10_h200_20261028.slurm).

This run tests whether the candidate geometry path recovers from the known
capacity failure and whether the sparse disk remains accurate through the
historical t=6–10 window. It is not a local-resolution comparison and does
not establish the GIZMO t=600 result.
