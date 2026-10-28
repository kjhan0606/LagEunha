# Matched native face-pressure force at t=6.5

The accepted aligned and 5% site-jitter bounded q2 snapshots at exact
`t=6.5` were copied with their original binary SHA-256
`3fc66a5f…` into an isolated A10 allocation, Slurm **437351**. The
[staging manifest](staging.json) records both parent snapshot hashes.
The binary's existing `GFS_GPU_FACE_DUMP_CALL=1` feature captured native
GPU face geometry, EOS inputs, and work for **the same 64 disk IDs** before
one short RK4 step to `t=6.5005`. The [selection](selected_ids.json) uses
eight seeded random IDs in each of six initial-radius bands plus 16 IDs
with the largest jitter native-minus-exact-orbit density residual. These
are selected sites, not a population-weighted sample of the disk.

The native calculations both exited zero; both [aligned](aligned/audit_short.json)
and [jitter](jitter/audit_short.json) 11-check short-replay audits pass. The
Slurm job itself is marked `FAILED 1:0` because its final Python
postprocessor lacked a helper import in the staged root. No GPU calculation
was restarted. The corrected auditor ran once on the saved native face
inputs and produced the [paired JSON result](force_probe.json). Its source
SHA-256 is recorded in that result. The copied raw TSV inputs and faces
retain their SHA-256 hashes in [copied_inputs.sha256.json](copied_inputs.sha256.json).

For each real directed face, the audit reconstructs the **active**
variational force

`F_i,face = (P_i-P_j) [A_ij/2 + L_ij(c_ij-(x_i+x_j)/2)/|x_j-x_i|]`

from the GPU's own face/input dump and divides the sum by fixed cell mass.
All 64 cells in each run have closed face area vectors. The maximum
native-work versus reconstructed `v_i·F_i-P_i Vdot_i,face` error is
`6.40e-25` aligned and `1.49e-23` jittered. Stored face pressure equals
the captured owner pressure exactly, and area vectors agree with polygon
endpoints within `9.13e-16`.

| Equal-ID group | Aligned median `|a_P|` | Jitter median `|a_P|` | Ratio |
|---|---:|---:|---:|
| 48 stratified random disk sites | `1.68150e-8` | `1.49516e-7` | `8.89` |
| 16 large-jitter-residual sites | `1.66207e-8` | `1.80078e-7` | `10.83` |

For the 48 random sites, the jittered radial and tangential pressure
acceleration RMS values are `7.657e-8` and `2.985e-7`; the aligned values
are `1.727e-8` and `8.991e-9`. Thus the jittered perturbation has a
prominent **tangential** component that can exchange orbital angular
momentum locally. The sampled median jitter pressure acceleration is only
of order `10^-7` in the run's units; this one-time force does not establish
the integrated orbital error or a physical viscosity.

An independent [pressure/geometry comparison](../KEPLER_Q2_BOUNDED_JITTER005_BALLISTIC_PRESSURE_T65_2026-10-27.json)
finds jittered native relative pressure fluctuation `2.831e-3` versus
`2.465e-3` predicted by each cell's adiabatic pressure on exact-orbit
Voronoi geometry, with mass-weighted Pearson correlation `0.981`. Together
with the near-zero signed entropy drift, this supports the hypothesis that
fixed-mass geometric area variation creates a pressure perturbation that
feeds back on the orbit. It remains an inference: a single force snapshot
does not distinguish source terms over `0<t<6.5` or prove a cell-spin cure.
The paired `t=10` runs remain the next evolution gate. Time-dependent
Laguerre volume control is a distinct candidate to test against this
measured force, subject to its weight/face conservation and runtime gates.
