# Guard against silently dropping extended checkpoint state

`kh_readdata` historically converted a mismatched stress-record stride by
copying only its common prefix. That drops the stress-resident thermal64 tag
and data and may cause a later thermal64 initialization to promote the
lower-precision mirror as though it were an untagged legacy state. A known
stress record converted to a non-stress target can likewise discard thermal64
or moment state.

The new reader rejects any stress snapshot whose stride differs from the
current build's stress-particle ABI. For a current-ABI stress snapshot read
into a non-stress target, it rejects a thermal64 or moment-tagged record
before discarding its stress portion. Legacy non-stress input and untagged
current-ABI stress-to-non-stress conversion retain their prior paths. The
behavioral tradeoff is explicit: older stress snapshots from another ABI now
require a deliberate version-aware converter; an unknown layout cannot be
validated merely from its byte count. Equal byte count is also not proof of
the same field layout, so a versioned checkpoint format remains needed.

The [guard test](test_kh_restart_guard_20261029.c) compiled and passed on
lageunha with the native `-DINTEL -DXYZDBL -DUSE_CUDA` type configuration:
legacy and matching stride are accepted; shorter and longer stress strides,
thermal tag and moment tag are detected. The modified `kh.c` also passed
`mpiicx -fsyntax-only` with the same defines. These are guard-level and
compilation checks, **not** an end-to-end native checkpoint restart or an MPI
test. The running Ada executable and frozen job scripts were not changed.
Before production use, build a separate native binary and exercise actual
legacy, tagged, incompatible, truncated and continuous-vs-restarted files,
including multi-rank behavior; the existing snapshot format has further I/O
validation limitations beyond this guard.

## Actual reader-path fixtures

The [reader harness](test_kh_restart_reader_20261029.c) links the actual
`kh_readdata` function and uses singleton MPI on lageunha. It stubs particle
migration and creates isolated one-particle fixture files in fresh `/tmp`
directories; it does not evolve the fluid or access an existing snapshot.
Both a CUDA-defined stress reader and a CPU non-stress reader accepted the
old header format, a new-format non-stress record, and an untagged current
stress record. Both rejected a deliberately incompatible stress stride with
exit 99. The non-stress reader rejected current-stride thermal64 and moment
tags with exit 99. The stress reader accepted those tagged records and the
test verified exact preservation of the thermal64 double value and moment tag;
moment reading used `GFS_MOMENT_INIT_ONLY=1`. These exercise the actual reader
branch and MPI abort in one process, not CUDA kernels, migration, real legacy
file compatibility across ABIs, truncated-file handling, multi-rank I/O or
continuous-versus-restarted evolution. Those production gates remain open.

## Truncated-file and header rejection

The native reader now checks every header, particle-record and required
time/timestep read; it validates positive particle count, legal stress flag,
minimum common-prefix stride and multiplication bounds. Before allocating
from the file's count, it checks that the file contains the claimed complete
records plus the required time, dt and previous-dt fields. The reader does
not demand an exact total file length, preserving legacy variations in later
footer fields. Truncation or an impossible header calls `MPI_Abort` instead
of continuing with partially read state.

The actual-reader harness was rebuilt on lageunha in CPU and `USE_CUDA`
configurations using singleton MPI. Old headers, current non-stress and
untagged current stress records remained accepted; thermal/moment tags were
preserved in the stress build and rejected in non-stress conversion. Both
builds rejected truncated headers, records, and required footer data, a
too-small prefix, a claimed two-billion-particle file with insufficient
bytes, an invalid stress flag, and incompatible stress stride. This checks
the CPU file reader compiled with the CUDA configuration; it does not run
GPU kernels, exercise corrupted-but-same-length data, validate snapshot
checksums, or prove a continuous-versus-restarted simulation. The existing
running binaries and frozen submissions were not modified.
