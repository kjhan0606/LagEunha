# Eight real coupled steps: serial versus MPI4

Both runs resumed preserved checkpoint 2305 at t=10.0020667693911, using
the previously verified experimental coupled executable (SHA256
fc2d1b5168dee1e531c44c2bfd77960f31f4e11e8f10680faa8b68679a25334e).
Separate new root: /gpfs/kjhan/LagForce/Codex_mpi_multistep_20261018.
The original trajectories were untouched. Runs used OMP2, CUDA nearest
neighbors, CPU force/HG, fixed masses, zero weights and no steering.

Both exited zero after eight accepted dt=0.0005 steps. Final snapshot time
is 10.006066769391104. All 36864 persistent IDs and explicitly compared
physical fields are bitwise identical, including position, velocity, mass,
density, pressure, internal energy, volume, weights, ceiling and sound speed.
All accepted states are finite with zero retries. Every RK stage mask is 15.
MPI4 records 48 reconstruction phases with consecutive epochs 1 through 48.

The final comparison and accepted ledgers are stored in
canonical_mpi_multistep_2026-10-18.json. Scripts:
run_mpi_multistep_probe.py and verify_mpi_multistep_probe.py.

This verifies repeated use of the coupled MPI path over eight actual steps.
It does not demonstrate migration crossings during this specific interval,
intermediate-stage rejection, active entropy-file output, CUDA hydro, dynamic
weights, or long-term physical accuracy. The physical interval is only .004;
this is a backend consistency check, not a replacement for t6-10/t20/long
Kepler tests. Existing single-rank physical instability remains unresolved.

Next: extend backend evidence to an interval with measured ownership changes
and repeated rejection, while assessing the original long-time trajectories
and the geometric/material-volume inconsistency independently.
