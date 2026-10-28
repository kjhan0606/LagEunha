# TCF full face algebra on independently reconstructed t6 geometry

For 11498 cells at r<1.9 about (3,3), SciPy Voronoi reconstructed all 68988
directed faces in the existing t6 dump. No selected face was missing or
infinite. Face areas differ from printed values by at most 4.915e-7 relative.
Using saved rounded Q, the independent assembly relative L2 discrepancies
against original cell totals are: force 7.134e-8, energy 6.371e-8,
volume rate 8.172e-13, directed heat 7.882e-8. These are consistent with the
limited precision of saved Q; this is not an exact production replay.

`Exam/tcf_full_face_probe.h` transcribes the default full TCF face assembly:
normal pressure, tangential pressure, dissipative force, energy work,
geometric volume rate and heat. On these reconstructed inputs,
`Docs/check_tcf_full_face_backends.cu` observed actual OMP teams 1 and 4;
CPU/OMP and CPU/CUDA outputs were all finite and bitwise identical.

This advances beyond the prior synthetic tangent and FACE_Z-only probes.
It still does not exercise original production code call sites, QBP gradient
construction, geometric kernels, per-cell parallel accumulation, MPI ghosts,
HG or RK4. The production GPU unsupported guard remains intact, and no
existing trajectory was restarted. No long-term physical pass is claimed.

Artifacts: canonical_tcf_full_face_backend_2026-10-18.json and
canonical_tcf_full_face_geometry_2026-10-18.json. Remote inputs, owners and
outputs: /gpfs/kjhan/LagForce/Codex_tcf_full_face_probe_20261018.

Next integration gate: recover the same per-cell totals from owner-indexed
face outputs under MPI partition/reduction and threaded accumulation, before
connecting this experimental kernel to production. Upstream QBP/geometry
parity and physical stability remain separate required gates.
