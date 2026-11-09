# Approved KH continuation, 2026-11-09

User explicitly approved existing KH and matching GIZMO continuation to t=10, intermediate density/conservation comparison, increased MPI ranks, and a larger simultaneous Slurm calculation. Wait for approval before any next scientific stage. Do not resume the broad autonomous goal or release previous held jobs.

Artifacts: `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/kh_t10_approved_20261109_v1`.

- 128×256: existing native common C=.2 trajectory resumes at1.35, unchanged executable and physics, one MPI rank and8OpenMP threads. GIZMO resumes the synchronized exported1.35 state with snapshot-restart flag2 and4MPI ranks×2threads. This rank change resets the GIZMO timestep hierarchy; it is not a binary-restart continuation. Both are running on LagEunha. Native is substantially slower.
- 256×512: same continuous KH initial condition, newly sampled131072cells, t0→10. Frozen portable native executable, same C=.2 and cell-scaled dt cap. Slurm457379 was submitted for H200 and then updated in place to h100, gres=gpu:H100:1, at user's request. One GPU,16CPUs,32GB,48h. Native remains1MPI×8OMP; GIZMO4MPI×2OMP. The source explicitly prohibits multi-rank native material evolution; do not remove that guard or claim that native MPI was implemented.
- Slurm is pendingQOSGrpGRES: part_h100 has aggregate GPU limit8, currently8allocated, although10physicalGPUs exist. User's normal-QOS per-user limit8 is not the blocker (user currently uses1GPU). Existing held jobs remain held.
- Compare density and conservation at2.1,4.7,9.2,10. Intermediate native restart segments are<=.39 to preserve the existing .4-window contract. Existing checkpoints are preserved. Physical-solver code is unchanged. A native H100 short portability gate precedes the high-resolution evolution.
- `worker.py` runs native C/CUDA and native GIZMO; Python is orchestration/postprocessing only. `monitor.py` observes and writes `n*/progress.json`, audits GIZMO snapshots, and creates matched-time comparison figures as both endpoints become available. No new physical model or subsequent test is automatically launched. Files show latest saved and observed times; do not describe queued n256 as running.
- GIZMO intermediate snapshot energy includes output-time prediction convention. The four-rank final restart exporter was verified bitwise against the archived one-rank synchronized export. At t10 it separately reports true synchronized endpoint conservation.
- Current native128 process1321174, worker1321002; GIZMO launcher1321155, worker1321003; parent1320975. Monitor PIDs are in progress.json, not durable truth; re-check process identity.

This is a running stage, not completed accuracy qualification. No current-method KH endpoint beyond1.35 existed before this approval. Paper McNally IC is different; do not call this direct reproduction of the published Figure16/17.
