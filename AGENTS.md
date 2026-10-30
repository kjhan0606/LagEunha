# Kepler project instructions

- Do not use the `/scratch` path for build files, temporary files, logs, or simulation outputs. Keep them in a designated `/gpfs/...` run directory, and inspect the submitted job script for `/scratch` references before submission.
- Run Kepler simulations through Syntax Slurm. Syntax is a shared login node; use it for editing, small compiles, submission, and result review, not for simulation runs or `mpirun`.
- Preserve existing run outputs. Use a fresh, job-specific build directory under the selected `/gpfs` run directory and refuse to overwrite files from a previous attempt.
