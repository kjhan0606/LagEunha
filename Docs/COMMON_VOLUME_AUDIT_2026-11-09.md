# Common finite-step volume audit, 2026-11-09

The previous goal turn integrated and qualified opt-in numerical solver reuse. This turn moves the existing immutable-ID stage audit from Cylinder-only ownership into a shared native include, samples the same realized volume rate used by pressure work in the periodic/reflected-wall RHS, and enables the audit for the already supported KH/RT/Sedov/Gresho paths. Noh remains excluded until its independent RHS is connected. This is diagnostic infrastructure for the common endpoint inverse, not another physical model or a new correction term.

The first RT diagnostic run stopped because its phase-3 gravitational-potential endpoint evaluation was mistaken for a third SSPRK2 stage. The revised sampler records phases 1 and 2 only and accepts the phase-3 diagnostic only after both samples exist. The failed private run is preserved; it was a diagnostic wiring error, not a demonstrated RT physics failure.

Qualification used existing saved states, four steps of dt=0.000125, without restarting the long calculations. KH: t=0.004 to 0.0045, 2048 cells. RT: t=3.001 to 3.0015, 8192 cells. Sedov: t=0.1 to 0.1005, 1024 cells. For each, all 13 checkpoint fields are bitwise equal between the original executable, new executable with audit off, and new executable with audit on. Each audit recorded eight stage samples and four endpoint comparisons. The existing Cylinder endpoint-inverse control is also bitwise unchanged. C/MPI source compiles and links to the existing CUDA operator; this does not qualify distributed MPI or 3D.

Measured maximum relative cell-volume discrepancy V_end - V_start - dt*(Vdot_1+Vdot_2)/2:

| Test | Maximum relative discrepancy | Maximum absolute summed stage volume rate |
| --- | ---: | ---: |
| kh | 6.36128759e-13 | 1.41922614e-16 |
| rt | 3.95351997e-07 | 2.19449296e-17 |
| sedov | 4.01096612e-06 | 1.79078801e-16 |

The largest absolute signed sum of endpoint discrepancies is below 1e-17. Global compatibility is therefore observed at roundoff in these intervals, while RT/Sedov have larger local discrepancies. This alone neither proves accumulated heating nor validates whole-test accuracy or morphology.

Next: connect the existing endpoint volume target to the periodic/closed-domain weight Jacobian, remove only its constant-weight gauge, retain the original full-row residual, and reject incompatible targets instead of subtracting their mean. Match this same volume update to pressure work and preserve fixed mass. Use the successful Kepler geometry/work consistency as the invariant; do not replace its tested baseline. Finish 2D before 3D. The GIZMO matched-accuracy total-time requirement remains unverified. The second paper is not being drafted before the first astronomy draft.

Evidence: `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_volume_audit_20261109_v1/results_v2/{comparison,summary,cylinder_parity}.json`; source preimages, private builds, commands, logs, checkpoint provenance, canonical hashes and patch are retained in the parent directory. The high-resolution fixed-wall Cylinder production binary/process was not changed.
