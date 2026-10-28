# Opt-in compact reconstruction implemented

Exam/exam.c adds GFS_QR_SMOOTH=1, default off. It multiplies each weighted
LSQ design/RHS row by (1-d^2/(6.25V))^2 inside the existing2.5sqrt(V)
radius, zero outside. This affects reconstructed velocity, pressure and
density, not just the velocity QBP input. Cell masses, Laguerre weights,
geometry force formulas, and particle motion laws were not changed.

The experimental mode requires cubic rank-revealing reconstruction.
Candidate overflow, invalid length scale, or insufficient rank aborts with
owner diagnostics; no silent wider-shell or quadratic fallback is used.
The actual QR operator export includes the taper in coefficients while
retaining physical source offsets. MPI/GPU-force coverage is not inferred.

Build: /gpfs/kjhan/LagForce/Codex_smooth_qr_build_20261018/eunha2.exe.
Separate exam.c/exam.o/library/executable preserve original running binaries.
Build commands inherit the previous CUDA-linked production build and change
only this source object; the production hydro calculation is still CPU.

## Terminal production validation

All three probes(off_plus,on_plus,on_minus) exited0 at preserved t3.4 states.
Default-off RHS differs from the older binary by at most2.71051e-20;
cross-build bitwise identity is NOT claimed.

For each on-branch,11498 current-r<1.9 cells were checked:

- Actual exported support matches complete point-cloud radius search:0 mismatches.
- All independently reconstructed cubic fits have rank9.
- Exported unlimited gradients match independent smooth fits within1.113e-12.
- Normalized cubic polynomial moment error<=2.687e-14.
- Original positions, velocities, masses and initial thermodynamic values match.

The first analyzer incorrectly compared unlimited operator export directly
with the limited stress gradients in the RHS file and stopped. Inspection
of alpha_vx/alpha_vy application corrected this comparison. Each branch has
22 limited velocity components; limiter collinearity residual<=5.652e-7,
consistent with rounded diagnostic fields. The large raw unlimited/limited
difference is recorded transparently in the JSON, not hidden or interpreted
as a failed polynomial reproduction check. No simulation rerun was needed.

Evidence: canonical_smooth_qr_probes_2026-10-18.json.
Scripts: run_smooth_qr_initial_probes.py, analyze_smooth_qr_initial_probes.py.

## Next actual-step gate

Four opt-in OMP4 runs from saved checkpoint743 at dt1e-4 and5e-5 were
submitted in /gpfs/kjhan/LagForce/Codex_smooth_qr_steps_20261018, session10037.
Runner check_coupled_generator_step.py now accepts --binary and --smooth-qr.
Its smooth results deliberately omit comparison with the legacy generator:
the generator must be reconstructed using the changed operator first.
No long-time, shock, MPI or CUDA-force stability claim follows.

Original coupled continuation to t4: reference,+,- each exited0 after760
steps since restart224. The five refined-to-t3p4 runs(session23450) remain
running. Original physical long runs remain separate from these diagnostics.

All four smooth one-step runs subsequently exited0, exactly one finite
accepted step each, no retries, unchanged masses. Inner departure-response
forward derivative is2.72878551e-7 at dt1e-4 and2.72806016e-7 at5e-5.
At the smaller step the legacy value is3.13787542e-7: the smooth candidate
reduces this local positive growth rate by about13.1%, but it remains
positive. This is NOT a stabilization result. A matching smooth-generator
budget and actual longer evolution are required next.
Evidence: canonical_smooth_qr_steps_2026-10-18.json.

## Matching smooth-generator comparison completed

Actual smooth hydro RHS and HG mode dumps were reanalyzed. Snapshot hashes
are identical to the states used for the existing pressure/kinematic terms,
so those terms can be reused. The HG generator is unchanged at these states;
the hydro remainder becomes more damping in the inner region.

Inner generator2.72739131e-7 versus forward derivatives2.72878551e-7 and
2.72806016e-7 gives relative errors+.05112% and+.02452%. Shear-band errors
are+.05268% and+.02529%; central errors-.09158% and-.04440%. These improve
approximately with the expected first-order forward-difference behavior.
Edge errors remain-.14666% and-.16543%; core signals~1e-15 have~1% relative
differences. Do not claim every region has converged or that this is RK4
order verification. Growth remains positive.
Evidence: canonical_smooth_qr_generator_2026-10-18.json and
canonical_smooth_qr_generator_comparison_2026-10-18.json.

## Physical comparison through t10 launched

Root /gpfs/kjhan/LagForce/Codex_smooth_qr_physical_to_t10_20261018.
Session82741, target10.002905917763403, MPI1/OMP4, dt cap.0025:

- legacy_continue PID3186276 resumes completed reference checkpoint984
  at t4.0029059 with the original executable; no t2-to-t4 legacy rerun.
- smooth_from_t2 PID3186283 starts the changed method from preserved
  reference checkpoint184 at t2.0029059 with the separate smooth executable.

Both use new directories. Dumps every2 time units plus final, estimated
snapshot storage345048060 bytes including source copies and allowance;
observed filesystem free286310656901120 bytes. Full environment and hashes
are in each case directory and root manifest. No initial-condition entropy,
mass, velocity or geometry reset is introduced at the method switch.
Adaptive stepping may produce slightly different dump/final times, which
must be recorded and used in physical phase-error calculations. These runs
cover the historical t6-10 failure interval, not the full long-time goal.
Retain original long runs and refined-response session23450.
Runner: run_smooth_qr_physical_to_t10.py.

## Physical evaluation pipeline and preserved baseline

analyze_smooth_qr_physical_snapshots.py checks complete snapshot sizes before
reading, validates IDs/masses/zero weights and finite positive states, then
caches actual-time physical metrics in root/physical_frames/case. It uses
the same original-IC material masks and np.linalg.norm radius calculation
as the prior t12 physical evaluation. This avoids mixing the different
roundoff-sensitive boundary memberships used by some response diagnostics.
It records wrapped phase, radial error, circular-velocity departure, exact
orbital-position error, relative pressure RMS and log-entropy RMS.

Preserved legacy t4.0029059 reference, inner background(804 cells):
phase RMS1.51835deg, relative-pressure RMS.338565, radial RMS.00359772,
circular-velocity RMS.00845082, exact-orbit position RMS.01035779.
Disk-core phase RMS is.00053755deg; relative-pressure RMS.00263528.
This is physical trajectory error, not a finite-perturbation growth norm.

The smooth method starts at t2.0029059 with already nonzero inherited inner
errors: phase.099168deg, relative pressure.024195, velocity.00041883.
These values describe its copied source state, NOT a smooth-evolution result.
No new physical-comparison evolution dumps existed at this analysis yet;
both live processes had begun accepting finite, zero-retry steps.
Evidence: canonical_legacy_t4_physical_2026-10-18.json,
canonical_smooth_start_t2_physical_2026-10-18.json.

## OpenMP and rejected-step recovery gates completed

Same smooth binary, t3.4 +/- states, accepted dt5e-5:

- OMP1 and OMP4 produce bitwise-equal id,x,y,vx,vy,mass,volume,den,
  pressure,ie,w2 in both branches.
- OMP4 forced final rejection of a dt1e-4 trial, followed by accepted dt5e-5,
  produces those same reader fields bitwise equal to the direct dt5e-5 run.
  Each rejection test has exactly one accepted finite step with retries=1,
  both exit0. Tests used new directories; no production trajectory rollback.

Full snapshot hashes differ; equality is scoped to the listed reader fields.
These are local new-option parity/recovery checks, not long MPI or CUDA
hydro validation. Per the established physical-first priority, no new MPI
expansion was launched while the single-rank long-time solution is unresolved.
Evidence: canonical_smooth_qr_threads_2026-10-18.json,
canonical_smooth_qr_retry_2026-10-18.json.

Live observation2026-10-17T21:36:52.873890+00:00: actual process handles
3186276/3186283 present. Legacy accepted54 steps to t4.1379059, smooth65
to t2.1654059, all accepted finite with0 retries. Restart-relative dE is
-3.90645e-9 and-4.48157e-9, respectively; these are not whole-trajectory
physical-accuracy measures or directly comparable same-time errors.
Evidence: canonical_smooth_qr_early_live_2026-10-18.json.

## Raw hydro conservation/thermodynamic checks

At preserved t3.4 +/- states, smooth TCF RHS has normalized summed force
residual<=1.24e-15 and summed total-energy-rate residual<=2.67e-17.
Interior identity dE-v.F=-P*dV+TdS closes to relative L1<=2.24e-13.
All221159 dumped directed faces per state have Q*du>=0 at printed precision.
These checks concern raw hydro before gravity, HG, or frozen-cell evolution
constraints; they do not establish accepted-step or long-time conservation.

The exploratory remote hydro_conservation.json incorrectly treated column10
dte as internal energy, double-counting kinetic work. Source inspection and
TCF derivation establish that dte is total energy. The original exploratory
artifact is retained but superseded by hydro_conservation_corrected.json;
its apparent energy residual is invalid. The corrected standalone script
is audit_smooth_qr_hydro_conservation.py. No solver change was made for this
analysis correction. Evidence: canonical_smooth_qr_hydro_conservation_2026-10-18.json.

## Same-time physical comparison procedure

After running analyze_smooth_qr_physical_snapshots.py on lageunha, run
compare_smooth_qr_physical_frames.py there. It uses cached complete-frame
audits, excludes the copied smooth source state, and compares only actual
times agreeing within1e-10. Duplicate legacy history/continuation frames
must have identical hashes. Region cell counts, fixed mass and zero weights
are checked before reporting each physical error and smooth/legacy ratio.
Comparisons are stored without overwriting under physical_comparisons.
Unmatched times remain pending; no interpolation or time substitution is
silently applied. The initial execution correctly excluded t2.0029059 and
produced no physical improvement claims because no evolved smooth frame
was yet available. This is analysis preparation, not a new physics result.

## Physical continuation stopped before first evolved dump

Actual smooth executable PID3186350 disappeared; its timeout parent3186283
was defunct. The log ends in MPI_Abort(2), not an observation timeout:
`[QR_SMOOTH] insufficient rank owner=539 rank=0 nuse=8 nc=140`.
Last accepted step724 has elapsed1.81, absolute time3.812905917763375,
finite=1/retries=0. Abort occurs during the next RK force evaluation4349.
Source inspection establishes that rank=0 here means QR was not attempted
because nuse8<ncol9; it is NOT a measured rank-zero design matrix.
There are140 searched candidates, but only8 inside the fixed2.5sqrt(V)
support. The strict no-widening/no-fallback guard intentionally stops.

This invalidates the fixed-support smooth option as a general continuation
method in its present form. It does not establish a physical blow-up or
that smoothing itself is ineffective. The new run directory contains only
the preserved t2 source snapshot; no evolved physical snapshot was saved
before this abort. Hence no same-time physical improvement claim is possible.
The runner waits for legacy first, so smooth exit.code can be absent until
legacy finishes; the actual process and MPI abort log establish termination.
Legacy executable3186332 remains live and must be preserved.

Next investigate a continuously varying reconstruction support with enough
geometric samples, including search coverage and conditioning, rather than
removing the guard or restoring abrupt shell switching. Neighbor-count-based
support is a candidate, not yet an approved mathematical closure: distance
order changes and rank degeneracy still need analysis. Preserve the failed
log and all original trajectories; no identical restart was launched.
Evidence: canonical_smooth_qr_rank_stop_2026-10-18.log.
