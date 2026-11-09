# Common vertex-positive primitive candidate: completed comparison and rejection

The preceding goal turn made progress by isolating the negative scalar reconstruction and producing a geometric positivity prototype. The user then identified an excessive cycle of local tests without decisive2D outcomes. This stage therefore implemented one unchanged common candidate, completed the existing1D suite and both predeclared native2D Sedov endpoints, and closed it on the full comparison. No additional shock-specific repair or subsequent candidate-flow run is introduced.

## Common implementation and execution

Artifacts are in `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_vertex_primitive_20261109_v1`. The new `gfs_material_positive_linear.h` scales a scalar gradient around the volume centroid, checking every supplied cell vertex. It is the same host/device function in the1D interval and2D polygon adapters, applied to BOTH density and pressure. There is no new physical floor or coefficient. A64-representable-step limit bounds only roundoff rechecking; failure is returned rather than tolerated. The generic dimension argument is not3D hydro validation.

The raw primitive pressure-offset QR, face-local scalar interval and tangent velocity pair are unchanged from the tested hypothesis. Exact Riemann, paired pressure force/work, mass and SSPRK2 remain unchanged. The helper refactoring reproduces all five prototype RHS arrays bitwise. The first1D build lacked the case header; that compile failure is preserved, the existing unchanged case header was copied, and no evolution used the incomplete build.

The native executable retains the separately qualified automatic search-capacity main object and the earlier inverse/accepted-step volume-norm agreement with unchanged1e-11 tolerance. Three archive objects were rebuilt together to retain the tangent-entry ABI: exam_gpu_extract.o, material_cuda_project.o, material_cuda_rhs.o. Original production sources/executables and the ongoing RT control were not replaced. Candidate executable SHA256: `97aa96d716e1784a5f8a11b0dd7fa7c7722a9b3b82f50fac33aa619ce3dc9e68`.

## Completed1D evidence

All nine original problems reached their original endpoints, with exactly the same IC text and case settings. Shu–Osher uses800 physical cells; the others use200. Every cell mass remains fixed, contacts remain at roundoff and the largest energy-balance error is4.04e-15. The cell-integrated reference norms are retained without substituting endpoint success for accuracy.

|Case|Original density L1 (%)|Candidate density L1 (%)|
|---|---:|---:|
|sod|0.56338|0.43159|
|blast|4.11098|4.30123|
|noh|0.40998|0.43754|
|shu|0.45933|0.30125|
|lax|0.48794|0.29576|
|doublefan|2.34847|2.24610|
|collision|0.34549|0.32346|

Blast peak density decreases9.0578→7.0251 toward the Athena reference maximum6.4417, while integrated density/pressure errors slightly worsen. Noh peak density instead increases4.9199→5.1345, and the wall-cell deficit remains. Double-rarefaction integrated pressure improves, but the central density error remains larger than the original control. These local limitations are not relabelled resolved.

The18 actual initial/final1D states embedded in weighted2D rectangular strips reproduce force, volume and energy rates with maximum scaled discrepancy4.56e-14. This checks interior dimensional reduction, excluding two outer cells with different stencils, not complete boundary or3D evolution. The actual native CUDA scalar+RHS objects agree with independent CPU evaluation to1.84e-15 overall; cached pressure discrepancy is2.64e-14. These implementation checks do not establish shock accuracy.

## Completed native2D outcome

Both64² layouts use exactly the preserved IC checkpoint, parameters, source radius, EOS, CFL and originalt=.1 target. Cell masses stay bitwise fixed, volumes and internal energies stay positive, and the actual geometric-volume criterion passes. Shared Ada CPU12,13, full native MPI-linked/CUDA executable, single rank. No repeated CPU gradient audit is included in timed evolution; the exact CUDA objects were independently verified before launch.

|Layout|Density L1 control→candidate (%)|Pressure L1 control→candidate (%)|rho=2 radius RMS control→candidate (%)|
|---|---:|---:|---:|
|glass|7.2095→11.2453|4.4504→6.6789|0.7635→0.9709|
|regular|10.9362→15.6575|9.5776→12.4219|0.9040→1.3348|

Both error norms and the measured contour roughness worsen in both layouts. Visually inspected unsmoothed actual-cell comparisons show larger low-density pockets and stronger shell irregularity in the candidate. Candidate maximum densities are5.60825(glass) and5.20097(regular), compared with4.32822 and4.79239 for the controls. Shared0..5 colors explicitly indicate over-range values. The analytic reference remains an ideal point explosion, while the matched ICs are finite-radius disks; this caveat applies to both methods.

Glass required1015 accepted steps/101 retries and188.998s; regular844/23 and129.064s. Old controls required844/16 and814/14, respectively (the authoritative control result JSON retains the exact count). The observed shared-GPU wall ratios are about.876 in each pair, but these are noncontemporaneous measurements at unequal accuracy. They establish neither a stable speedup nor the <=2xGIZMO goal. Total-energy relative changes are−1.45e-12 and−6.17e-13; good conservation did not prevent degraded morphology.

## Decision

**Reject this candidate as a common-method replacement.** Its improved admissibility and local continuity do not compensate for the native2D accuracy regression. Keep the original successful method and all negative evidence. Do not launch additional shear/Gresho/KH/RT or3D calculations with this candidate and do not add a per-test limiter, sensor, floor or coefficient to rescue it. `decision.json`, the two endpoint results and the density figures are the authoritative closure.

The full objective is still incomplete: accurate dimensionally common hydro, the remaining2D validation and references, later3D, and matched-accuracy <=2xGIZMO runtime are not established. Independent original RT continues. Further research must prioritize a measured evolved physical defect and the complete common-method outcome; do not treat another passing local reconstruction check as qualification. The first astronomy manuscript precedes drafting the second technical paper. No manuscript or remote repository was modified in this stage.
