# Native Q2 pressure-sensitivity experiment

This experiment tests whether pressure changes the short-time growth of the
dominant inner-ring angular-momentum mode on a fixed existing perturbation.
It does not replace the original GIZMO-like Kepler IC or prove long-time
stability. The preserved Q2 `t=30` checkpoint (55,351 cells) is the common
source; the baseline copy is byte-identical. The second copy multiplies
internal energy and pressure by four, sound speed by two, and consistently
scales stored entropy and pressure-gradient fields. This deliberately creates
a **different physical state**, not an ordinary restart of the original
solution. Positions, velocities, masses, other record bytes and metadata
remain unchanged. The [preflight auditor](audit_kepler_pressure_scaled_ic_20261029.py)
checked exact W0, mass, byte-level changes and EOS consistency; its maximum
EOS relative error is 1.01e-7, consistent with the stored float fields.

Both runs use the same historical native Q2 binary
`87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968`, same
parameters, fixed polygon, eta=.06, no steering and Re=infinity. The binary
contains an sm_90 image. Each independent Slurm array element requests one
H200, two CPU cores, 32 GB and two hours; each advances from t30 to
32.501708696164926 and writes one final snapshot, map and logs. There is no
dependency between array elements. A prewritten [comparison auditor](audit_kepler_q2_pressure_compare_20261029.py)
measures the same original inner-ring `0.2<=m/N<0.3` mode, energy relative to
each run's own initial state, mass and W0 invariants, positivity and entropy.
Additional native logs, geometry and wall checks are required before accepting
the endpoints. A pressure change also alters the background and sound-speed
CFL, so any result is a controlled sensitivity, not proof of a universal
`c_s/d` law.

The private staged directory is
`/gpfs/kjhan/LagForce/Kepler_Q2_PressureCompare_20261029`; the
[stage script](stage_kepler_q2_pressure_compare_20261029.py),
[snapshot writer](kepler_scale_thermal_checkpoint_20261029.c), and
[Slurm script](kepler_q2_pressure_compare_20261029.slurm) preserve the source
and verify frozen input hashes. Initial submission 440182 was canceled while
reported pending to make `SLURM_EXPORT_ENV=ALL` explicit for the MPI step;
the original staged script and hash manifest are retained with a
`before_export_fix` suffix. Track the replacement array and inspect actual
job status and outputs; do not infer success from submission. The H200
`sbatch --test-only` estimate was very late and is a scheduler estimate,
not a guaranteed start time.

Replacement array **440185** was submitted with `--export=NONE` and the
explicit MPI-step export; both elements were observed `PENDING (QOSGrpGRES)`.
The canceled array was confirmed `CANCELLED` at zero elapsed time. The
comparison auditor was exercised with a deliberately artificial mock in which
both cases point to the old unscaled t32.5017 endpoint; that checks file and
calculation plumbing only, and the x4 mock numbers are **not** a simulation
result. No new evolved results are available at submission.

An additional [quality auditor](audit_kepler_q2_pressure_quality_20261029.py)
checks the pressure-comparison receipt against independent reconstructed
power-cell areas for 50,514 original `r<2` cells, full-domain volume, density,
EOS, fixed polygon and W0 pressure modes, halo clipping, logged injection,
boundary/geometry failures, and each case's own energy reference. It reads
the log only through the saved endpoint's exact step marker: later failures
must not be attributed to an earlier checkpoint. Its reference-path mock
used the historical t32.5017 checkpoint and historical log for both cases;
all listed geometric and log checks passed, with active-cell area error
5.84e-8. This mock does not validate an x4 evolution or prove all production
acceptance criteria. The first mock attempt correctly exposed a JSON boolean
serialization problem, and another exposed the need to truncate later log
events; both were fixed and rerun before use.
