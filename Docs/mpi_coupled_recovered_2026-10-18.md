# Recovered full MPI coupled-step evidence

Existing output directory Codex_mpi_coupled_step_20261017 already contained
completed serial/np1/np2/np4/retry2/retry4 outputs. No simulation was relaunched.
Its verify.py was inspected and run on lageunha. All five comparisons succeeded.

From t=10.0020667693911 to t=10.0025667693911, all 36864 persistent IDs match.
The comparator verifies time/dt, finiteness, and explicitly listed fields:
x, y, mass, vx, vy, den, pressure, ie, volume, w2, w2ceil, csound (as present
in the saved layout). Every reported field is bitwise identical to serial.
MPI1/2/4 each records six reconstruction phases and RK stage mask 15. Forced
final-trial rejection on MPI2/4 records twelve phases, one retry, and one
accepted step. Resulting physical fields match the direct serial smaller step.
Maximum compared accepted-ledger difference is 3.02e-13 from reductions.

Configuration: OMP_NUM_THREADS=2, OMP_DYNAMIC=FALSE, shared geometry,
FACE_Z/QBP, frozen-consistent split HG, GFS_MPI_COUPLED_EXPERIMENTAL=1.
CUDA nearest-neighbor path is enabled; TCF force and HG remain CPU.
Saved and current executable SHA256 both equal:
fc2d1b5168dee1e531c44c2bfd77960f31f4e11e8f10680faa8b68679a25334e.

This supersedes claims that coupled MPI RK/retry has never been exercised.
It does not prove repeated-step integration, all intermediate rejection
locations, active entropy-file output, dynamic weights, GPU TCF, or physical
long-term stability. It tests the experimental entry path, not unrestricted
MPI support for every configuration. Diagnostic flags/guards remain intact.

Full recovered results: canonical_mpi_coupled_step_recovered_2026-10-18.json.
Authoritative original outputs and verifier:
/gpfs/kjhan/LagForce/Codex_mpi_coupled_step_20261017.

Next MPI gate: repeat actual coupled steps from an existing disturbed-state
checkpoint, comparing serial and distributed states at the same accepted
times, ownership crossings and fresh phase stamps. Preserve original long
trajectories. Avoid repeating already completed static geometry/RHS probes.
