# Native U64 frozen-W0 Kepler build

The optional `Makefile.native-u64` compiles the original CUDA, MPI and OpenMP
program with explicit double internal energy, thermal RHS and cell volume.
The precision state is an appended particle-owned payload. The ordinary
build retains its576-byte particle record. The optional build uses696bytes,
version2 checkpoints and the unchanged576-byte prefix. Weight position,
momentum, persistent inertia and RK increments are represented explicitly.
The qualified long mode requires exactly frozen zero weights and momentum.

Use a FRESH separate checkout and a Slurm compute allocation. Existing
archives are refused and library compilation is forced, preventing a mixture
of576-byte and696-byte objects. No source or object from an older experiment
is needed. Build prerequisites are Intel oneAPI2025.3, Intel MPI2021.17,
GNU13, CUDA13.0.2 and the existing FFTW3f installation at
`/home/kjhan/local`. The CUDA object currently targets H200/Hopper sm90.

```bash
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export BUILD_JOBS=2
make -f Makefile.native-u64
```

Request at least2CPUs and16GiB. Compilation needs no GPU. Set TMPDIR, TMP and
TEMP to a directory under/gpfs. The executable is
`.native-u64/eunha2_native_u64.exe`. The source-build script compiles all25
local support archives, all10 native solver objects and the main program.
The compile-time layout check verifies both the default576-byte record and
the explicit696-byte record, including the576/608/672 appended offsets.
`Exam/NativeU64` contains the qualified native source and exact header
versions previously selected by the experimental include path. The root
`eunha.h` adds only conditional payload declarations, leaving default layout
and ordinary solver sources unchanged.

Long mode requires `GFS_NATIVE_U64_W0_EVOLVE=1`, a version2/features7 accepted
checkpoint, MPI1 and the existing fixed polygon contract. The reader,
RK stages and writer refuse nonzero W/pi or a nonfrozen payload. The mode
never resets a weighted checkpoint. `GFS_WEIGHT_NATIVE_RK4_VALIDATE=1`
continues to permit only the separately bounded weighted short experiment.
The expensive acceptance audit can print every100steps through
`GFS_U64_AUDIT_EVERY=100`. Positivity, energy/work, rollback/cache and wall
ledgers are still checked on EVERY step. Accepted cumulative thermal-store
roundoff is printed at normal completion with a since-restart scope.

Checkpoint format conversion does not recover precision lost in an old
trajectory. Start a thermal comparison from the unchanged original IC with
exact float-to-double promotion, preserving all masses, positions, velocities
and original thermal values. Resume accepted checkpoints and carry the entire
initial log lineage. Do not use the ordinary legacy reader for696-byte data.
The new native codec and mode guards are mandatory.

## Evidence and limits

The staged fresh source build used jobs446036 and446037. The first compiled
all libraries and native units, then stopped on a missing FFT include in the
main compile. The second reused only those newly compiled local objects and
compiled/linked the main program with the corrected includes. No simulation
was rerun to repair a build or reporting failure.

The saved parity package is
`/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_fullsource_parity_20261101_v1`.
All51checks pass against the earlier qualified native binary with actual
OMP8 and MPI1. The comparison covers20accepted steps, an actual forced
full-payload rollback, and a fresh2-step continuation. All19registered
physical fields,120explicit payload bytes, footer, accepted time/wall
records and cumulative roundoff match bitwise. Process-local pointers make
whole-file hashes differ. The default/native layout checker separately
prints576/696bytes. Native comparison outputs pass the strict C codec.

The already-running independent physical job446033 retains its original
qualified executable. No running job was restarted or swapped to obtain this
engineering result. The source-build parity does not establish long stability,
weighted long-run accuracy, MPI multi-rank physics, other hydrodynamic tests
or AREPO-relative performance. The original t55 reference has28passing
aggregate physics checks but fails3of9wall qualification checks. All such
failures remain authoritative.
