# Native q2 GPU tessellation resource screen

The accepted `t=10→20` aligned W0 A10 log reports a `density_copy_A` time of roughly `150–200 ms` for one GPU tessellation call, versus single-digit milliseconds for the other device substeps. Four tessellation calls occur per RK4 step. Kernel A therefore remains the main measured target for a long `t=600` calculation. The [completed same-checkpoint Ada comparison](KEPLER_Q2_GPU_TESS_ADA_MATCHED_GATE_2026-10-27.md) found the 3×3 seed slower; the pending H200 duplicate 437439 was canceled without running.

Read-only `cuobjdump --dump-resource-usage` of the isolated v3 and v4 CUDA binaries reports the following **sm_90** resources for `voronoi_tessellate_kernel`:

| Binary | Candidate seed | Registers/thread | Stack/thread | Shared memory |
|---|---:|---:|---:|---:|
| v3 `7358c549…` | 5×5 | 128 | 3,280 bytes | 0 |
| v4 `945864b0…` | 3×3 + certified expansion | 128 | 3,296 bytes | 0 |

The sm_80 version of v4 also reports 128 registers and a 3,280-byte stack per thread. The source declares both an input and a temporary polygon, each with 48 double-precision vertices and support IDs, as well as other per-thread data. The reported stack is consistent with substantial thread-local storage, although this resource report alone does **not** measure spills or prove that local-memory traffic dominates timing. Seed3 changes candidate enumeration but does not reduce the stack allocation. Its speed and exact geometry have now been tested at one accepted checkpoint on Ada.

The first in-place clipping experiment failed native geometry at 151 cells, while a limiter-mode specialization with the original clipper gave only a 0.4–1.8% Kernel A gain despite a 752-byte stack reduction. These actual tests show that stack size alone is a weak acceleration target. Any next candidate must retain the 48-corner abort and support-plane certification and be compared against an accepted checkpoint in both W0 and weighted modes; candidate/clip work deserves direct measurement.
