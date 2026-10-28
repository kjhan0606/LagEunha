# Candidate: use the previous stage's faces as an early clipping seed

The [seed-bin exclusion experiment](KEPLER_WEIGHTED_SEED_BIN_CERTIFICATE_2026-10-27.md)
was bitwise equivalent in short W0 and weighted tests but approximately
2% slower on the tested A10. It left the expensive initial scan for
directional anchors intact. A different candidate is to use real neighbors
from the immediately preceding tessellation as early clipping seeds, then
still execute the complete current seed and adaptive searches. The previous
faces are a **performance hint**, never the only neighbor search.

At each current stage, a cached neighbor plane must be built from the
neighbor's **current** coordinates and weight. Each such plane is then an
ordinary current constraint, even if that neighbor no longer owns a face.
The complete current search and geometric checks remain necessary to find
new neighbors and validate the final result. Changing clipping order can
change roundoff and transient polygon complexity, so algebraic inclusion
alone does not prove bitwise equality or avoid the 64-corner abort.

The current `FaceCSR` stores neighbor indices, not the previous owner and
neighbor IDs. Raw indices cannot safely survive reordering, migration or
halo changes. A prototype would need to capture the previous owner ID and
neighbor IDs before the particle SoA is overwritten, and check them against
the current owner and neighbor indices before using any cached plane. If
identity, bounds, finite positions or weights fail, use the current full
directional-preclip path. Duplicate periodic images can share an ID; until
image identity is also certified, disable this shortcut for any padded or
multi-rank case. The no-padding W0 and weighted 2D controls can test a first
candidate, but they do not satisfy the final MPI requirement.

Retain the existing 32-plane near-wall safeguard. For an interior cell, a
valid cached seed must also close the initial box; if it does not, discard
the partial polygon and execute the existing directional preclip. A rare
intermediate corner overflow or any geometry failure must be recoverable by
the original per-cell path, not accepted as a truncated face. Rebuilding
GPU buffers or a changed particle count must invalidate the cache. Stream
ordering must ensure previous faces and IDs are captured before the current
upload and before face buffers are overwritten. These requirements are not
implemented yet.

CPU Slurm **437541** measures independently reconstructed neighbor-set
persistence over already completed one-step W0 and weighted replays. That
observation can tell whether a warm start is worth implementing; it cannot
verify CUDA buffer lifetime, cache validity, MPI, conservation, or speed.
Any implementation still needs matched physics, face geometry, bitwise or
quantified error, fallback and timing tests at early, distorted, weighted,
boundary and multi-rank states before promotion.

Job **437541** completed `0:0`. In the [W0 jitter t=6.5→6.501](KEPLER_NEIGHBOR_CACHE_STABILITY_2026-10-27/w0.json)
replay, 99.920% of the final directed disk-neighbor relations were already
present at the start; 99.070% of disk cells retained their entire neighbor
set. In the [weighted jitter t=0.825→0.826](KEPLER_NEIGHBOR_CACHE_STABILITY_2026-10-27/weighted.json)
replay the corresponding figures are 99.938% and 99.264%. For W0, 3.05%
of cells in the inner initial-radius band `0.5–0.75` changed neighbor set,
versus 0.46% in `1.5–2`; the two bands contain different cell counts and
physical scales. These two isolated steps make warm starts worth exploring,
but do not prove a causal link to heating or predict cache hit rates at other
times, resolutions, MPI decompositions or longer steps.

CPU job **437542** sampled 1,024 disk cells in each of the same two replay
intervals. Its [geometric prototype](audit_kepler_cached_seed_geometry_20261027.py)
clips the current local box first with either the existing eight directional
anchors or the previous snapshot's real neighbors evaluated at current
coordinates and weights. In the [W0](KEPLER_CACHED_SEED_GEOMETRY_2026-10-27/w0.json)
and [weighted](KEPLER_CACHED_SEED_GEOMETRY_2026-10-27/weighted.json) samples,
all cached seeds closed the local box and used 5.99 and 5.98 planes per cell
on average, compared with eight directional anchors. The median cached to
directional squared circumradius ratios were 0.9913 and 0.9935. The mean
number of seed-bin candidates inside that initial radius was about
`8.87→8.03` for both cases; these counts are upper bounds before the
subsequent streaming cuts.

This sample supports a feasibility investigation, not a speed result: it
does not execute or time the CUDA kernel, test all cells, exercise the native
intermediate-corner fallback, or validate MPI and halo images. The substantial
potential benefit would be avoiding the first full directional-candidate scan,
not merely skipping empty bins, but its actual cost and safety remain to be
measured.

## Active-kernel safety audit — 2026-10-28

Read-only audit of the Q1 running executable's archived source overlay (`q1/dteta_no_cap_from_t1p5_to65_20261028/source_overlay/Exam/exam_gpu.cu`) sharpened the required fallback:

- In raw-geometry mode, the starting polygon is a local square with `halfbox=min(boxsize/2,cellsize)`. After the complete search, the kernel treats any negative `poly_neigh` support or any coordinate still on `halfbox` (within `32*DBL_EPSILON`) as an artificial-box touch. It then doubles the box and repeats, up to the existing hard failure. A cached warm start can skip the directional prepass only after its clipped polygon has no artificial support under this same predicate; otherwise it must restore the initial square and run the exact existing prepass.
- `gpu_clip_voronoi_halfplane(..., exact_intersections=1)` is not a recoverable trial primitive: its raw append path increments the global geometry-error counter and returns `-1` on corner-capacity failure, and the orchestrator later aborts on that counter. The helper does commit its temporary polygon only after the full clip succeeds, so a candidate-specific exact clip helper can safely return a local failure without mutating the polygon or poisoning the global error counter. On such failure, restore the initial square and use the existing path.
- Therefore a first candidate needs a recoverable exact-clip helper plus explicit per-row cache validation (owner and neighbor IDs, real faces only, in-range SoA indices). Preserve the complete current 5x5 seed scan, adaptive expansion, final geometry checks, and box-doubling logic. A failed check, ghost/periodic image, MPI padding, changed particle count, or incomplete cached closure must use the baseline path.
- Capture the previous `pindx` array after the current scatter on the same CUDA stream, without overwriting the previous ID snapshot before Kernel A has validated and consumed it. Do not enable the shortcut for MPI/periodic-image cases until image identity is certified.

This is an implementation constraint, not a speed or physics result. The active W=0 t=1.5-to-6.5 job continues with its archived binary and source hashes unchanged. The rejected seed3 and seed-bin-pruning candidates remain excluded.
