# Actual compact injection and original-mass restriction completed

2026-10-22. Implementation under the selected
`ASTRA_COMPACT_SPACE_AFTER_COST_FAILURE_20261022.md`. The preceding goal turn
made progress by completing the paired-output audit and structural preflight;
this turn makes progress by constructing reusable native coefficient and mass
data. Neither turn proves the t6–10 heating/stability objective. No physical
time evolution was performed here. The preserved old accepted step5 remains
at t=0.00011263377801805954.

## Native injection, not merely a topology count

`build_compact_injection_20261022.cpp` constructs the actual constrained
polynomial injection and two dynamic initial-velocity residual columns.
It reads the saved canonical reference, source ownership and conditioned
native basis. It preserves the fixed core r<0.5 choice and shared trace
constraints. Its four new rank shards bind the old fine identities, sparse
polynomial rows, binary80 residual columns, compact local/global identities,
owners and initial coefficients, with complete payload FNV checksums.

Results in `material_compact_injection_v1_mpi4`:

| Item | Actual result |
| --- | ---: |
| MPI ranks | 4 |
| Compact global scalar DOF | 1,581,183 |
| Protected source cells | 5,524 |
| Coarse/fine interfaces | 240 |
| Maximum polynomial injection row width | 9 |
| Maximum initial position coefficient difference | 1.1102230246251565e-16 |
| Maximum initial velocity coefficient difference after reconstruction | 0 |
| Actual native preparation time | 7.09525226 s |
| Output size, shards plus result | 882,273,610 bytes |

| Rank | Unique owned compact DOF | Ghost compact DOF |
| --- | ---: | ---: |
| 0 | 395,910 | 0 |
| 1 | 395,295 | 615 |
| 2 | 395,295 | 615 |
| 3 | 394,683 | 1,227 |

Native MPI ownership counts sum to the actual compact dimension. The initial
x/v owner/ghost values were pulled and compared bitwise. The position
difference from rounded reference interpolation was reported, not erased by
changing the old geometry or generating thermal context. Zero initial
velocity coefficient difference does not prove all future field accuracy.

`material_compact_injection_io.h` implements full checksum/ABI/identity reads
and distinct transpose operations for local element loads versus uniquely
owned assembled fine loads. Ghost-replicated assembled fine loads must not
be transposed over all rows, because that would count shared rows twice.

## Real binary80 mass congruence and retained weak RHS

`assemble_compact_mass_20261022.cpp` read those complete injection shards and
the authoritative stored binary80 M. It assembled **T^T M_old T**, combining
the four children outside the protected core and retaining protected fine
elements. This reuses the original mass integration; it performs no new
quadrature and does not create an independent thermal or spin reservoir.

Results in `material_compact_mass_v1_mpi4`:

| Item | Actual result |
| --- | ---: |
| Compact element blocks | 351,112 |
| Coarse parent blocks | 262,440 |
| Protected fine blocks | 88,672 |
| Maximum local polynomial support | 13 |
| K_old under stored M | 8.78884830216854831195 |
| K_new under congruence | 8.78884830216854831021 |
| Signed K difference / original U | -2.0077827583e-15 |
| Signed L difference / (2K) | 0 |
| Signed retained-load work difference / U | -2.8904601003e-19 |
| Signed retained-load torque difference / (2K) | -8.8476807855e-23 |
| Old-measure velocity difference squared | 0 |
| Old-measure position difference squared | 5.6507787718e-37 |
| Actual native assembly time | 4.68226314 s |

The raw residual Gram has diagonal about2.27974501008419355e-15 and
off-diagonal about1.47213270124e-32. Suggested column norms are about
4.77466753825e-8. These norms permit future border normalization without
altering the represented physical field; normalization must also change the
corresponding coefficients consistently. Positive residual Gram alone does
not establish positivity of the full Schur complement.

The saved weak RHS is **T^T(M_old a_saved)**, with a_saved read from the
original converged initial acceleration. It inherits that original solve's
true residual (about9.8e-17). It is not a bitwise saved exact force vector,
an average of old accelerations, or a fresh force query at the slightly
different rounded compact initial map. This distinction remains explicit
before an actual compact initial operator is evaluated.

Stored matrices are symmetric by evaluating each upper-triangle contraction
and storing the same result in its reflected entry. No diagonal clipping,
mass lumping or post-hoc kinetic-energy correction was used. K/L/work/torque
comparisons check this assembled result against the original stored M.

## Solver implication and next implementation

At hanging interfaces a protected fine triangle can reference thirteen
compact polynomial masters, although its polynomial functions span at most
the ten-dimensional P3 space. Its local polynomial block can therefore be
singular despite positive diagonal entries. The original per-element SPD
Cholesky preconditioner cannot simply be reused with this changed support.

The next operator must assemble the global polynomial mass action and handle
the two global residual columns as a rank2 border/Schur system. It must
verify the actual Schur positivity and original1e-16 true equation residual,
then evaluate the unchanged full initial pressure-curvature/thermal/work
gates. The existing CUDA/DD vector and MPI halo work remains available for
this smaller system. These preparation times are not RK4 throughput and do
not satisfy the t1/t10 performance gate by themselves.

## Integrity and execution scope

Both native commands exited0; their retained exec sessions55428 and73479
completed and were collected. The source/input/binary integrity records
before and after each run all report OK. The new mass completion record
checks all four actual file headers and binds SHA256 of the new outputs.
The mass assembler itself checked every injection payload checksum.
Each wrapper and its environment copy were made read-only before launch;
neither running source nor wrapper was edited. No new scheduler job, new
agent, optional guard experiment or independent physical trajectory started.

The implementations use the saved native basis and MPI/OpenMP. Python was
used only for manifests/result bookkeeping, not physical evolution. The
original source repository, old accepted checkpoints, generating entropy,
shared global context and all earlier evidence remain preserved.

Pending: initial curvature/pressure and full true-residual solution; actual
compact CUDA force and cubature cost; continuous/global map validity and
material shear lifetime; real t1/6/8/10 density/pressure/orbit/rotation/energy
behavior; W0 dynamic Voronoi and weight/mass contracts; shocks, boundary
matching and the longest correctly normalized GIZMO comparison.
