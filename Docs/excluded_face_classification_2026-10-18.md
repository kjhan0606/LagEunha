# Excluded faces concentrate at the zero-weight density interface

Classify every previously excluded unique internal face, using saved t2
generator states and the specified Wdot=target probe.

| Weight fraction | Excluded | Both endpoint W and Wdot exactly zero | Overlapping donor intervals | Undercovered intervals |
|---|---:|---:|---:|---:|
| .25 | 945 | 945 | 921 | 24 |
| 1 | 955 | 955 | 935 | 20 |

Midpoint radii range approximately .498--.609, around the inner disk edge.
For overlaps the common situation is two donor intervals each covering the
same entire face. The face coincides with an original donor boundary and
cannot be interpreted as a single continuous donor trace. Undercoverage is a
separate clipping/degeneracy issue and must not be described as proven physical
two-state overlap.

Median absolute normal velocity jump is .03447/.03456, maximum .088653.
Maximum |Pj-Pi|/(Pi+Pj) is .01643. These jumps are diagnostic states; no shock
classification follows from them. The largest normal jump includes cells
841--957 near r=.5154.

## Consequence for the candidate

The earlier retained-face calculations did not cover many of the zero-weight
faces at the important density interface. Therefore good conservation on the
retained subset must not be presented as coverage of the most relevant physical
interface. Positive density fits or compact segment algebra do not fill that gap.

Local W_i=W_j=Wdot_i=Wdot_j=0 is NOT by itself the user's global, persistent
equal-relative-weight reduction condition: other generators have nonzero
weights/rates, and the diagnostic examines one instant. Nevertheless these
faces are a concrete transition/interface case any full candidate must handle.
Do not force a zero flux here by an unsupported local switch, nor choose a
donor arbitrarily. Such shortcuts can destroy legitimate transport or break
the approved reduction rather than prove it.

Next build both one-sided traces on the coincident donor boundary and separate
the exact contact case (equal pressure and normal velocity) from general
velocity jumps. A physical interface solver and the original fixed-mass branch
must then have a consistent mathematical connection. Resolve the20--24
undercovered faces as a distinct numerical geometry gate.

No new physical solver or trajectory is authorized by this diagnostic alone;
existing user authorization remains in force, but scientific acceptance is
unmet. The original runs remain live: at2026-10-17 18:30:58UTC baseline t15.706640,
fixed-HG t15.509351, from-t0 time-cap t11.745497. No restart or cancellation.

Evidence: `canonical_excluded_face_classification_2026-10-18.json`.
Script: `audit_excluded_face_classification.py`.
Remote: checkpoint-remap directory, `excluded_face_classification.json`.
