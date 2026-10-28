# Experimental native W=0 vertex certificate

The files in `Exam/` are byte-preserving overlays on the frozen native clip
sources recorded in `MANIFEST.json`. No production executable was replaced.
The CPU-only Slurm build 439062 compiled the three affected translation
units and linked all four CUDA targets (SM80/86/89/90) in 86.86 s. Its
candidate SHA-256 is
`a22e0e8ebea0eca807ab3762d95bf037184b7784218a42fd2eb1ae7a24526bfc`.

`GFS_GPU_VERTEX_HALO_CERT=1` can select the new certificate only for one MPI
rank, the Kepler IC, kappa=0, exact zero weights and the standard periodic
padding callback. Other configurations retain the old test. Both force-stage
tessellation and accepted-state volume refresh receive the actual padding
rectangle; its upper bounds are not inferred from a rounded-up grid size.
The candidate-search and half-plane clipping algorithms, stress/gradient
arithmetic, flux/work, RK4 and fixed masses are unchanged in source. Their
numerical parity after recompilation still requires the executable qualifier.

`gfs_vertex_halo.h` compares an upward bound on squared vertex radius with
a downward bound on squared clearance from the actual rectangle. Exact
tangency and uncertain margins are rejected. CUDA uses explicit directed
rounding operations described in [NVIDIA's double-precision intrinsic
reference](https://docs.nvidia.com/cuda/archive/12.8.2/cuda-math-api/cuda_math_api/group__CUDA__MATH__INTRINSIC__DOUBLE.html).
The host test uses outward `nextafter` bounds. This bounds arithmetic on
the represented polygon; it does not prove correctness of the preceding
floating-point polygon construction. The former native completeness and
capacity checks remain, and independent geometry/force parity is mandatory.

The direct diagnostic `GFS_GEOM_FORCE_DUMP_MAXCALL=4` is read-only and off
by default. It saves the first four RK4 force calls' complete face CSR and
per-cell acceleration, thermal-work rate, entropy-work rate and timestep.
The comparator normalizes neighbor IDs and their periodic image displacement
because halo width changes padding indices. It still requires exact vertex
coordinates, face ordering, supports/images and force/work arrays. It rejects
incomplete files and wrong images even when generator IDs agree.

The C++ boundary tests cover exact tangency and invalid inputs plus 100000
near-boundary examples. All 50308 accepted examples satisfy independently
evaluated long-double squared inequalities. Four dump-format tests and ten
common gate tests pass, including rejection of force changes, wrong periodic
images and changes in owned/face counts. These tests do not establish CUDA
runtime parity or long-time physical stability.

## Submitted native experiment

GPU job 439065 has five cases: a canonical/new-executable qualifier with the
option off, then a control at halo40, an enabled halo40 control and an enabled
halo20 candidate. All start from the accepted t=2.6549183918680868 state;
the three full cases stop exactly at t=2.66 and export all four first-step
force calls. Exact endpoint fields, accepted-step trace, native physics gates
and every energy ledger must pass. The smaller halo must actually exercise
the new certificate. Wall measurements include the initial read-only dumps;
production speed without dumps remains a later check. About 1.2 GB of output
is budgeted. No early-time trajectory is restarted.

The pending host-padding/MPI/local-sync job 439029 is independent: it retains
its previously frozen executable and auditor. The new vertex experiment
keeps both padding reuse and local alias refresh off. The updated common
auditor permits differing padding counts only when explicitly requested by
the halo-width comparison; existing comparisons still require exact counts.

The running long trajectory 438957 remains on its original halo40 executable.
Neither t>6 nor GIZMO's t=120/t=600 density-field endpoints, the longer prose
orbit claim, conservative cell-spin evolution or weighted production parity
has been established by this work.

The build directory on GPFS contains all 67 compiler-declared local
dependencies, libraries and commands. It preserves the original library as
`libexam.input.a`; `rebuild.slurm` restores that input before repeating the
original preflight. The initial submitted script, logs and binary remain
unchanged. No build was repeated just to preserve this provenance.
