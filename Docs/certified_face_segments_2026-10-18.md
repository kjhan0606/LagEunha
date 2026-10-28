# Complete selected-face segmentation with certified one-sided donors

Replace the ambiguous polygon-line intersection partition by the lower envelope
of original donor squared-distance functions along each straight weighted face.
After removing the common u^2 term, all functions are affine in face coordinate u.
Find their crossing coordinates, select nearest donors on each open interval,
and determine left/right limits from normal directional derivatives for ties.

Candidate search is bounded geometrically: the face owner supplies a maximum
distance bound along the segment; all potential nearest donors lie within that
bound plus half face length of the midpoint. A KD-tree supplies this finite set.
Both interval endpoints are checked against every candidate. Affine distance
differences then certify the entire interval within the stated floating tolerance.
Adjacent intervals with the same donor pair are merged.

| Weight fraction | Selected unique faces covered | Segments | Discontinuous two-sided segments | Maximum segments per face |
|---|---:|---:|---:|---:|
| .25 | 3705/3705 | 5371 | 1239 | 3 |
| 1 | 3705/3705 | 5365 | 1233 | 3 |

No donor lies outside the r<.7 patch. Partition length sums are one to printed
precision. Maximum nearest-owner certificate residual is6.93e-17; moving-face
kinematic identity residual<=4.35e-17. The previous945/955 exclusions are no
longer needed for this selected-face geometry. Endpoint portions with third
donors are represented explicitly, not overwritten by the face owners.

This certifies the selected internal faces used by the diagnostic, not all
faces/boundaries of the global simulation. Floating tie tolerances are not an
exact-arithmetic proof and require broader degeneracy tests before production.
The observed maximum3 is not a universal memory bound under repeated evolution.

## Next physical interface gate

Stored records now provide face endpoints, common normal speed coefficients,
segment bounds, and both donor IDs. This allows exact contact and uniform-state
controls without arbitrary donor selection. General jumps still need a
physically consistent interface flux including pressure work and total energy;
the original fixed-mass reduction cannot be established from geometry alone.
No physical solver, production checkpoint, MPI ownership or time stepping was
changed. Long Kepler stability remains unproved.

Full face data are stored remotely as `certified_face_segments.json` in the
checkpoint-remap directory, with path and SHA256 in the repository summary.
Evidence: `canonical_certified_face_segments_2026-10-18.json`.
Script: `audit_certified_face_segments.py`.
