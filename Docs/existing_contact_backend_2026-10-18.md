# Existing CPU/CUDA contact-state helpers agree on exact contacts

Inspected `Exam/exam.c`, `Exam/exam_gpu.cu`, and the documented 3D copy
`Exam/gfs_riemann.h`. The HLLC helper returns star pressure and normal velocity,
not a complete Eulerian/ALE mass/momentum/energy flux. Do not treat it as a
ready-made replacement for the new transport model or the original fixed-mass
force path. The 3D header was read, not executed in this test.

Extracted the existing CPU hll_star_state, hllc_face_2d and rest-frame wrapper,
and their CUDA counterparts, stripping comments but preserving arithmetic.
Use postype=double, phase1=0, original pressure-ratio constants. Source SHA256
values and probe generation are recorded. No production source changed.

4096 manufactured exact contacts cover differing densities including .01/1,
pressures1e-5 or1, common normal velocities approximately +/-26.6 and differing
face-relative velocities. Compile on lageunha with GNU13/CUDA13, FMA and host
contraction disabled; launch an actual CUDA kernel.

- CPU one-thread versus requested OpenMP four-thread outputs: bitwise equal.
- CPU versus CUDA: bitwise equal, maximum difference0.
- Star-pressure reference error<=1.70e-21.
- Common normal-velocity reference error<=1.39e-17.
- All outputs finite; executable exits0.

OpenMP was configured with omp_set_num_threads(4); actual team size was not
instrumented in this first probe. This limitation does not affect the reported
output comparison but must be fixed before a strong parallel-execution claim.

## Scope and next step

These checks support reuse of the existing contact-state arithmetic. They do
not cover general unequal pressure/normal velocity, fallback/rarefaction
guards, phase1 extreme faces, full production kernels, MPI, or full ALE flux
assembly. The previously derived exact contact transport still requires
upwind density/tangent selection, pressure work and correct frame transforms.
The helper alone does not supply these quantities.

Next instrument actual OpenMP teams and compare general face-state outputs
on the certified original donor pairs, with branch coverage and frame/side
transform checks. Do not change the active fixed-mass solver merely because
the contact subset agrees. Long physical Kepler stability remains open.

Evidence: `canonical_existing_contact_backend_2026-10-18.json` and
`existing_contact_backend_provenance_2026-10-18.json`.
Generator: `build_contact_backend_probe.py`.
Remote probe/build/result: `/gpfs/kjhan/LagForce/Codex_existing_contact_backend_20261018`.
