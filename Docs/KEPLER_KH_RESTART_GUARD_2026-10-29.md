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
