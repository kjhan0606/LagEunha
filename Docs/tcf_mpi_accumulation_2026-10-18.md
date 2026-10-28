# MPI/OpenMP accumulation of reconstructed TCF face contributions

Executed on lageunha with gnu13/13.2.0 and openmpi4/4.1.5. The diagnostic
uses the previous 68988 reconstructed t6 face inputs and transcribed full
TCF kernel. MPI_Init_thread requested and received at least FUNNELED.

All twelve combinations passed: MPI ranks 1/2/4, actual OpenMP teams 1/4,
partition by owner ID or by face index. Global face counts are exactly 68988.
Thread-private cell buffers are reduced locally, then MPI_Allreduce sums
rank-local cell arrays. All outputs are finite. Largest absolute discrepancy
against sequential CPU cell totals is 2.0816681711721685e-17 (reference maximum
0.0059238251734637613). Prior CUDA face outputs summed sequentially produce
the same reference. Differences reflect changed summation order.

This is a real MPI collective test, but every rank reads the full saved input.
It does not test production partitioning, ghost construction/exchange,
geometry, QBP gradients, multi-GPU execution or per-rank CUDA device assignment.
No production MPI/CUDA correctness or physical stability pass is inferred.

Code: Docs/check_tcf_mpi_accumulation.cpp.
Results: canonical_tcf_mpi_accumulation_2026-10-18.json.
Remote binary: /gpfs/kjhan/LagForce/Codex_tcf_full_face_probe_20261018/mpi_probe.
Compiler flags: -O2 -fopenmp -ffp-contract=off. Each run used
mpirun --bind-to none -np N ./mpi_probe T with N=1,2,4 and T=1,4.

Remaining integration gate: inspect and compare production per-rank input
states, neighbors and gradients before considering activation of the TCF GPU
path. Existing unsupported guard and running trajectories remain unchanged.
