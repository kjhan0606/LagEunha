# Remaining native transfer execution

Goal: original CUDA/MPI/OpenMP Kepler stability past t6--10 and ultimately the
GIZMO longest run. This CPU binary80 reference acquisition is a required
transfer stage, not an evolved solution or GPU acceleration result.

Use the fixed465 remaining overlaps,84489 panels,21629184 positive Gauss16x16
queries and the same source thermal context, density/eps/profile and original
all12 tau. Original68526 and the completed four whole targets are excluded.
No new atlas, arrangement, event inverse or implicit real branch solves.

Geometry is the original reference chart. Cached regular/oblique values
retain full uncertain root/free intervals. The direct native input passes
both binary80 lambda centers and the original chart lambda enclosures to
cache(), alongside original Cartesian and weight balls. Gradients and all
qnode/edge caches are unchanged. Legacy callers keep their default path.
No lambda normalization, force steering, heat reset or extra energy reservoir.

Before physical queries, freeze oracle sites: first Gauss point for every
slot/type/original chart (whole,radial,angular,oblique), plus one saved failed
legacy input per source/kind/edge mask. These same stream points are compared
against the original Q128 same-function field17 and all12 moment enclosures;
their physical queries are appended once. Pointwise native forward-error,
edge/tail/J and heat positivity guards remain mandatory for every query.

Output bound:21802217472 positive bytes plus84489*24640=2081808960 commit
bytes, approximately23.89GB, plus small oracle/summary files. GPFS free222TiB
was checked before launch. Outputs are unique, existing runs never replaced.
OMP16/MPI1 runs manually on lageunha; login syntax performs only editing and
bounded reads. CPU native walltime is measured, not attributed to GPU speed.

Each source overlap has canonical panel/iu/iv order. Flush positive data
before appending a per-panel commit with exact next-query index, rolling
FNV64 byte checksum and full64-level binary80 reduction tree. Controlled
failure also saves the partial-panel prefix and tree. Every numerical point
is stored once. FNV is an accidental-corruption check; final SHA256 manifests
bind files and source/evaluator/input policies. Files are fflush-ordered;
power-loss durability beyond normal filesystem completion is not asserted.

This first-acquisition driver does not yet implement resume. Any continuation
must validate frozen inputs/arithmetic-policy and positive prefix against
the saved commit/tree/checksum, then resume exactly the next query. Never
infer success or restart only from file size, lockfiles or an observation
timeout. Successful targets or point prefixes are not recomputed.

Each overlap uses the same value/error packet for future source outgoing and
target incoming. Current execution emits readonly target candidates only;
full384 source closure and full24 target transactional state commit are
separate required steps. Failure of an overlap prevents its target from
publishing. Actual native/error/tree, saved uniform analytic remainder and
final R80 conversion must satisfy original target tau and physical m,U
positivity before any target candidate passes. No caller-state commit here.

No accepted RK4 step, evolved heat/restart support or t>6 stability is claimed
from this run. Those remain the next implementation stages after closure.
