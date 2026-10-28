# Bounded Kepler weight-rate solver: device-local iteration screen

The accepted native near32 source was reconstructed in the isolated Git
worktree `/tmp/Kepler_volume_wdot_gpu_20261027` on `syntax` from base commit
`472e7fcc3419b060d8fa0c74f0961f6c1f9f188b` plus the 42-file native
and 7-file bounded overlays. Every copied file passed the archive SHA-256
check. The [reconstruction script](reconstruct_kepler_volume_gpu_probe_20261027.py)
refuses a second overlay and records the accepted binary hash
`3fc66a5f...`. This worktree is isolated from the live CUDA runs and is
the intended location for a weighted GPU diagnostic.

## Jacobi versus the existing 16-step PCG reference

The [read-only auditor](audit_kepler_q2_bounded_jacobi_rate_t65_20261027.py)
uses the actual accepted `t=6.5` power diagrams, static target-volume
weights, and exact face-area Jacobian. It solves the same
`J Wdot = -Vdot_fixed` equation as the earlier CPU reference. The Jacobi
iteration is `Wdot <- Wdot + omega * (rhs-J Wdot)/diag(J)` and needs
only per-cell face rows, with no global dot-product reduction. The
[aligned](KEPLER_Q2_BOUNDED_ALIGNED_JACOBI_RATE_T65_2026-10-27.json) and
[jittered](KEPLER_Q2_BOUNDED_JITTER_JACOBI_RATE_T65_2026-10-27.json)
receipts compare four relaxation factors and 8/16/32/64 iterations.

| `dt=0.00193983`, inner-band p99 area error | Aligned | 5% site jitter |
|---|---:|---:|
| Frozen weights, full Euler power diagram | `1.299e-4` | `2.013e-4` |
| Jacobi16, omega 0.8, **linearized** | `7.157e-6` | `7.216e-6` |
| Jacobi16, selected best omega, **full Euler** | `3.056e-5` (omega 1.0) | `3.156e-5` (omega 0.8) |
| PCG16, full Euler | `3.098e-5` | `3.103e-5` |
| PCG16, **linearized** | `2.856e-6` | `2.606e-6` |

The 16-step Jacobi linear residual is larger than PCG16, but the
**rebuilt one-step areas** reach the same `~3e-5` floor in both ICs.
This means the cheaper, synchronization-free iteration merits a native
GPU implementation check. It does **not** establish convergence of the
linear solve, a safe long-time update, or a GPU speedup. Omega was
screened on these same states; use one fixed value, likely 0.8, in an
independent restart before treating it as a production parameter.

## Independent device implementation check

An [isolated CUDA kernel](kepler_gpu_jacobi_probe_20261027.cu) now builds
each active owner's fixed-weight area rate and the conductances from a
directed power `FaceCSR`, then performs 16 device-local Jacobi passes at
`omega=0.8`. It leaves particle weights and the hydro force unchanged.
The [CPU exporter](export_kepler_q2_weight_rate_csr_t65_20261027.py)
provided the actual accepted `t=6.5` aligned and jittered power-face
graphs: 49,038 active cells and 297,114 directed faces per case. A10
Slurm job **437363** completed `0:0`. Its [verified receipt](KEPLER_Q2_GPU_JACOBI_VALIDATION_T65_2026-10-27.json)
records a synthetic two-cell/Dirichlet/unsupported-wall check and both
full-graph comparisons:

| CUDA versus independent CPU Jacobi16 | Aligned | 5% site jitter |
|---|---:|---:|
| Maximum absolute `Wdot` difference | `1.82e-15` | `1.58e-15` |
| Relative vector L2 difference | `1.19e-10` | `1.12e-10` |
| One device call, including temporary allocation/copy | `6.21 ms` | `5.48 ms` |

The [opt-in attachment](apply_kepler_gpu_jacobi_diagnostic_20261027.py)
was compiled in a separate source worktree and linked into a complete native
executable. It reads a persistent ID mask, calls the kernel once on the
**native tessellator's** device `FaceCSR`, and writes `Wdot` by stable
particle ID without changing the force or RK4 state. A10 Slurm job
**437365** completed `0:0` on aligned and jittered accepted `t=6.5` restarts.
The [native receipt](KEPLER_Q2_NATIVE_GPU_WDOT_PROBE_T6501_2026-10-27.json)
verifies that all ten recorded physical fields at `t=6.501` are **bitwise
equal** to the original-binary replay in both cases. All inactive rates are
exactly zero. Native GPU versus an independent CPU power-geometry Jacobi16
reference differs by at most `1.98e-15` aligned and `1.55e-15` jittered;
relative vector L2 differences are `1.28e-10` and `1.21e-10`. The probe
call, including scratch allocation and host copy, took `1.964 ms` and
`1.947 ms` on the actual native face graphs (333,111/333,086 directed rows).
This validates the **read-only stage attachment and rate arithmetic**; no
weight update or weighted hydrodynamic evolution has yet been performed.

## Native attachment constraints

The existing GPU tessellation leaves directed `FaceCSR` and `ParticleSoA`
on device before `getAccVoro2DBlend_GPU` launches its force kernel.
Each real owner's directed face provides face length, owner-to-neighbor
distance, centroid, generator position and velocity. Those suffice to
assemble `Vdot_fixed`, the diagonal `sum L/(2d)`, and one matrix-vector
product `sum (L/(2d))*(rate_i-rate_j)` per Jacobi iteration. Core and
exterior rates must remain zero, and periodic/MPI padding must map rates
from their real owner or be explicitly treated as fixed boundary values.
Do not classify active disk cells by their **current** radius inside the
generic solver; persist a mask/target state across the RK4 and restart
chain.

Before evolving weighted gas, ensure the same stage `Wdot` drives the
face-normal geometric speed and thermal work; derive the weighted
pressure force from the same power geometry. The current `GFS_W2_LINEAR`
path has the correct normal-speed term but is restricted to positive
hydro-kappa weights, while `GFS_DYN_W2` is CPU-only and solves a different
face-speed objective. Neither can be enabled unchanged for this test.
The subsequent [pressure-weight power audit](KEPLER_Q2_T10_AND_WEIGHT_POWER_GATE_2026-10-27.md)
establishes an additional condition: a bare fixed-mass update leaves
`dE/dt=-P^T J Wdot`. A global work-neutral rate changes the candidate
Jacobi16 rates by only 0.72%/0.32% in the aligned/jittered one-step
screens, but its local physical law and native finite-step conservation
are unverified.
The W=0 fixed-mass Voronoi reduction, pairwise conservation, positivity,
and per-stage timing against the original CUDA binary are release gates.
The true end target remains the published long-time GIZMO disk, rather
than the present closed-wall polar control.
