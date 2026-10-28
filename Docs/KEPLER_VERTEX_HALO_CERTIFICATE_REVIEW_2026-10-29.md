# W=0 vertex-disk halo certificate: review and first reference result

This is a geometry/computational-cost investigation. It changes neither the
ongoing Kepler trajectory nor its physics. Astra reviewed the mathematics in
one focused consultation; the driver implemented and executed the audit.

For a generator p and a bounded convex candidate Voronoi polygon with
vertices v, a missing generator q cannot cut the polygon if
`|v-q|^2 - |v-p|^2 >= 0` at every vertex: that difference is affine in v.
Consequently, if every disk centered at v with radius `|v-p|` lies inside
the rectangle D covered by the available generators, no generator outside D
can cut the polygon. The rectangle enclosing the vertex disks has bounds
`min_v(v_x-|v-p|)`, `max_v(v_x+|v-p|)` and the analogous y bounds.
Each disk lies inside the old owner-centered radius-2R disk, so the new
test is sufficient and can accept cases rejected by the old bound.

Astra identified prerequisites that remain mandatory before any native change:

- D must be the rectangle actually covered by complete halo data, including
  MPI corner images. An allocated linked-cell grid rounded upward is not
  evidence of generator coverage beyond the requested padding.
- The adaptive search must have processed all potential cutters inside D;
  the external certificate does not repair an incomplete internal stencil.
- The polygon must be finite, convex and free of active artificial-box or
  capacity truncation. Periodic generator images and vertices must use the
  same unwrapped coordinates.
- A tie is harmless for area in exact arithmetic but can affect face identity.
  Negative margins must not be accepted by adding an epsilon. The CPU audit
  uses a positive diagnostic margin, not a formal floating-point error proof.
- Restrict this proposal to W=0. Keep the existing clamped Laguerre guard for
  nonzero weights until a separate sufficient condition is derived.

## Read-only accepted-state reference

CPU Slurm job 439042 completed in 23 seconds, without a GPU. It read the
accepted t=2.6549183918680868 checkpoint (SHA-256
`c34a846dedaf2e01a8e5548de2aa091f98db53f761424ce6cd931107ff0fc4a4`)
and reconstructed all 147456 cells using nine periodic tiles. The periodic
image range passes the numerical vertex-disk coverage test with a minimum
margin of 5.0287. Reference areas sum to exactly 36 at the reported precision;
maximum relative difference from stored native float volumes is 5.9582e-8.
Peak process RSS is about 2.05 GiB. No trajectory state was advanced.

| Halo cells | Old bound: unproven cells | Vertex bound: unproven cells | Minimum vertex margin |
| --- | ---: | ---: | ---: |
| 12 | 280 | 128 | -0.22129 |
| 16 | 232 | 0 | 0.02871 |
| 20 | 180 | 0 | 0.27871 |
| 28 | 0 | 0 | 0.77871 |
| 40 | 0 | 0 | 1.52871 |

At this accepted state even the old halo28 test passes, with a small minimum
margin of 3.7439e-5. Therefore this result **does not show that the four
failures at the later geometry evaluation were false positives**. It shows
that a substantially smaller halo can be sufficient for this reference
geometry under the vertex criterion. The failed predictor geometry itself
has not been reconstructed here.

The next native experiment should compare the unchanged halo40 path with an
opt-in W=0 vertex certificate using the actual padded rectangle and a smaller
halo (20 is a candidate with more margin than 16). First require an old/new
flag-off executable qualifier; then verify effective edge geometry and
neighbor/image identity, per-cell force/work and accepted-state parity over
the former failure interval. Compare wall time on the same allocation. Do
not promote based on area totals or this CPU reference alone. The current
halo40 long run and the pending padding-reuse/local-sync comparison stay
intact. This work establishes neither t>6 stability nor intrinsic-spin or
dynamic-weight evolution.

Artifacts: [stage receipt](KEPLER_VERTEX_HALO_CERTIFICATE_STAGE_2026-10-29.json),
[reference result](KEPLER_VERTEX_HALO_CERTIFICATE_RESULT_2026-10-29.json),
[audit source](audit_kepler_vertex_halo_certificate_20261029.py), and
[three geometric acceptance/rejection tests](test_kepler_vertex_halo_certificate_20261029.py).
