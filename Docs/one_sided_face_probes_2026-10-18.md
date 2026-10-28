# Two-sided donor probes and endpoint certification

For every excluded face, query the original zero-weight donor partition at
u=.25,.5,.75, offset to each side by epsilon times generator separation.
Use epsilon=1e-5,1e-7,1e-9 and compare the last two scales.
All945/955 faces give stable donors matching the two weighted face owners at
these three interior samples, including the24/20 previously undercovered cases.
This alone is NOT proof of whole-face coverage.

Add a stronger endpoint test: compare owner distance squared with the nearest
original donor distance squared at both endpoints. Squared-distance differences
between any two generators are affine along a straight segment, so endpoint
nonnegativity certifies nearest-owner status throughout, to floating tolerance,
when the two owners share the original bisector.

| Fraction | Excluded faces | Interior probe success | Whole-face endpoint certificate |
|---|---:|---:|---:|
| .25 | 945 | 945 | 913 |
| 1 | 955 | 955 | 922 |

The32/33 remaining faces have portions near endpoints belonging to other
original donor cells. Maximum squared-distance gap is1.93e-5/5.93e-5,
far above the roughly5.7e-14 scale of the endpoint tolerance. Interior sampling
would miss these portions. Do not assign a single pair of traces to the entire
face based on midpoint tests.

## Next geometry/interface step

Split each ambiguous face at all relevant original donor-bisector crossings.
On each open subsegment determine the left/right limiting donors; a subsegment
may be single-donor or coincide with a donor discontinuity. Endpoints have
zero measure and need a consistent convention, not an invented physical state.
Use analytic distance-difference bounds to certify intervals rather than relying
on a fixed small offset everywhere. Then evaluate physical two-sided interface
states. Sampling is useful diagnosis, not the final flux selection algorithm.

This classifies geometry only; it does not choose a Riemann/contact solver or
prove fixed-mass reduction. No physical evolution, production checkpoint,
CUDA/OpenMP/MPI path or long-time stability was changed or validated.

The initial endpoint run failed JSON serialization of a NumPy boolean-derived
integer before writing results. Explicit bool conversion fixed reporting; no
numerical criterion changed. Failed script remains remote as v2; completed
script v3 writes `one_sided_face_probes_v2.json`. Earlier sampling-only output
is retained as `one_sided_face_probes.json`.

Evidence: `canonical_one_sided_face_probes_2026-10-18.json`.
Script: `audit_one_sided_face_probes.py`.
Remote: checkpoint-remap directory.
