# Separate accepted TCF and HG entropy diagnostics

Exam/exam.c now supports GFS_ENTROPY_SPLIT=1 with GFS_TCF_SACC, restricted
to MPI1 with split HG and TCF. Output retains ID and existing total first,
adding TCF RK quadrature and HG split log1p increments. Each component is
saved/restored alongside total entropy on rejected trials. Existing output
format is retained when the flag is off. MPI coupled path still rejects SACC.
These accumulators start at process initialization, not at the original IC
time when restarting; prior entropy history is not reconstructed.

Built separately at /gpfs/kjhan/LagForce/Codex_entropy_split_20261018, using
the prior build commands and unchanged supporting libraries. Tested one t10
step at dt=.0005 in three cases: component flag off, on, and forced final
rejection of dt=.001 followed by dt=.0005 acceptance. All listed physical
fields match the same-binary off case bitwise. Existing total entropy matches
on/off bitwise; both separated components match retry/direct bitwise. Maximum
total minus components is 6.74e-19. The off case still enables legacy SACC;
this gate specifically tests the new component flag's passivity.

An initial strict comparison against the previous executable FAILED bitwise
equality: 2 vy values differ up to 5.55e-17, 229 pressures up to 1.08e-19,
241 internal energies up to 1.06e-22, 121 sound speeds up to 2.78e-17. All
other listed fields match. This cross-build difference is preserved in
canonical_entropy_split_old_binary_difference_2026-10-18.json; its exact
compiler/source cause is not established. Do not claim cross-build bit parity.

New-component verification: canonical_entropy_split_2026-10-18.json.
Verifier: verify_entropy_split_probe.py. No live trajectory was altered.
No t0-to-t12 entropy split or long-term physical accuracy is established.
Next: use a preserved evolved checkpoint for a matched multi-step component
measurement, recording its restart origin explicitly. Verify the source gate
and compiler change separately before any production binary replacement.
