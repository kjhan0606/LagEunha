# Active FACE_Z dissipative algebra: isolated backend comparison

The t6 recorded face inputs supply 221159 directed rows. The standalone
`Exam/tcf_face_z_probe.h` evaluates FACE_Z with QLAW=1, using recorded fi/fj,
qin and cell velocity difference. It does not recompute QBP gradients.

Actual OpenMP team sizes were 1 and 4; outputs agree bitwise. CUDA is not
bitwise identical: maximum absolute difference is 1.734723475976807e-18.
An independent NumPy expression agrees to 3.469446951953614e-18 maximum.
No directed dissipative heat is negative.

Against the original rounded printed Zq and Q, relative L2 differences are
2.3617304114185652e-7 and 2.1139274656712694e-7 respectively. This is consistent
with recorded input/output precision, not a machine-precision production replay.

This establishes consistency of this isolated algebra on sampled inputs only.
It does not validate full pressure force/work, geometry, reconstruction, HG,
MPI accumulation, RK4 stages, or physical long-term stability. The active
production TCF path remains CPU-only; its GPU unsupported guard is unchanged.

Artifacts: `canonical_active_face_z_backend_2026-10-18.json`,
`canonical_active_face_z_analysis_2026-10-18.json`; remote binaries and inputs:
`/gpfs/kjhan/LagForce/Codex_active_face_z_backend_20261018`.

Live observation at 2026-10-17 18:52:52 UTC confirmed original PIDs and working
directories: baseline 2286374 at t=15.7599147650, fixed-HG 2295268 at
t=15.5422148073, from-t0 dt cap 2466162 at t=11.8990036456. All last accepted
states were finite with zero retries. No restart or cancellation was performed.
These facts do not establish physical accuracy of the trajectories.

Next backend gate must cover the actual full TCF shared-face force and energy
assembly, including orientation reversal, before attempting production GPU
integration. Scalar FACE_Z parity alone cannot justify enabling that path.
