# Native bounded q2 Kepler at exact t=10: aligned passes, 5% jitter fails

The original fixed-mass, `W=0` CUDA+MPI+OpenMP A10 continuation used the
accepted exact `t=6.5` checkpoints, with no restart from the original IC.
Both native solvers reached exact `t=10`. The aligned job **437328** finished
`COMPLETED 0:0` in 54:24 and passed all 17 physics checks. The jitter job
**437344** reached exact `t=10`, but its auditor exited `1:0` because the
`disk_ballistic_residual` check failed. The threshold is `1e-3`; the measured
mass-weighted absolute disk residual is `1.196684e-3`. Its other 16 checks
passed. This failure is retained in the
[jitter native receipt](jitter_failed_audit_t10.json).

| Metric at t=10 | Aligned | 5% site jitter |
|---|---:|---:|
| Mean absolute disk density change | `1.96891e-4` | `2.37446e-3` |
| Exact-orbit tessellation prediction | `1.94290e-4` | `1.47192e-3` |
| Native minus exact-orbit density residual | `6.03158e-6` | `1.19668e-3` (**fails**) |
| Mean absolute disk log-entropy change | `1.19396e-5` | `2.03302e-5` |
| Signed disk log-entropy change | `-8.46237e-6` | `-5.07674e-7` |
| RMS site-position error against exact orbit | `7.05295e-7` | `2.12454e-5` |
| Total-energy relative error | `8.10306e-9` | `8.07655e-9` |

The jitter residual rose from `3.44253e-4` at `t=6.5` to `1.19668e-3` at
`t=10`, while its exact-orbit geometry prediction stayed near `1.47e-3`.
The aligned residual remained `6.03e-6`. Every initial radial disk band
has a larger jitter residual: `1.59e-3` at `0.5≤r₀<0.75`, declining to
`8.37e-4` at `1.75≤r₀<2`. The inner jitter/ aligned ratio is 55.0;
the whole-disk ratio is 198.4. The small signed entropy and total-energy
errors show that this gate failure is not evidence of large **net** thermal
heating. The increased site-position error and dense alternating pattern
are consistent with a geometry/pressure feedback, but these endpoints
alone do not prove its force mechanism.

The [paired density field](density_exploratory_t10.png) uses common color
scales. Its top row displays surface density, and its bottom row displays
`Σ/Σ₀−1` in the disk. It rasterizes each saved cell density at the nearest
generator; the numerical metrics come directly from complete native
snapshots and independently reconstructed exact-orbit cell areas. The
[paired exploratory receipt](pair_exploratory_t10.json) and
[aligned](aligned_radial_tails.json)/[jitter](jitter_radial_tails.json)
radial receipts retain the failed-check label and parent snapshot SHA-256.
CPU-only job **437366** produced these diagnostics and completed `0:0` in
11 seconds. Its predecessor **437362**, which required both native jobs to
pass, was canceled while still pending with `DependencyNeverSatisfied`.

The aligned success releases the separately audited `t=10→20` job
**437343**. Its byte-verified accepted parent checkpoint SHA-256 is
`3714b691ba996f9725db07df2383b4a806f18a2264b8a875b05ce1b6ac97f9a3`.
The jitter state is for diagnosis and must not be promoted as a successful
long-time validation. This bounded polar-ring IC and wall are distinct from
the paper-matched GIZMO long endpoint.
