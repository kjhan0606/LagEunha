# Saved remap donor trace versus affine reconstruction

Interpretation qualification: see `remap_small_weight_limit_2026-10-18.md`.
Exact donor trace preservation is not a necessary criterion for all convergent
affine reconstructions. This audit alone cannot reject that entire class.

Actual t2.0029059178 conditional remap, not a new physical simulation.
Script: `audit_remap_face_trace.py`.
Evidence: `canonical_remap_face_trace_2026-10-18.json`.

The saved states contain aggregate moments but not donor fragments. Reconstructed
old and weighted polygons from the same checkpoint and saved weights, and used
the saved density-fit coefficients, without refitting or changing initial data.
Each new edge was clipped against old donor polygons. Ambiguous coincident
donor boundaries or incomplete interval coverage were excluded explicitly.

On retained edges, integrate two velocity traces with identical donor density:
(a) original piecewise constant donor velocity, (b) receiver affine velocity
recovered from its stored moments. Use the receiver mean velocity as a fixed
diagnostic face speed and measure mass, relative momentum and relative kinetic
energy flux in that frame. This is not the actual evolving Laguerre face speed.

| Weight fraction | Retained / eligible one-sided faces | Inner mass flux relative L2 difference | Inner relative K flux relative L2 difference |
|---|---:|---:|---:|
| .25 | 5448 / 7165 | .8504 | .9788 |
| 1 | 5439 / 7165 | .7132 | .8827 |

Inner relative momentum component differences are .945--.959 at fraction .25
and .811--.859 at fraction 1. Edge relative kinetic flux differences are near
one. These are ratios of facewise L2 norms, not a global energy loss fraction.
No cancellation into a closed-cell balance is claimed. Faces are one-sided
records and may include both sides of a shared physical edge.

Order8 versus order12 line integration differs by at most 6.64e-15 absolute;
sampled density remains above .00941. This checks quadrature consistency, not
the full clipping implementation or physical accuracy of either trace.

## Consequence

The residual internal velocity information is not negligible for face
transport in this diagnostic. A recovered affine slope alone is not an exact
representation of the remapped field, even though its volume moments satisfy
the previous positivity and energy split checks. Do not integrate the affine
headers as a complete closure or silently thermalize the residual.

The donor trace itself is not a validated physical reference solution: it is
the represented conditional lifted field. This comparison establishes loss
of representation, not a Kepler instability mechanism or a preferred solver.
There is no Riemann interface solution, pressure work, gravity, time update,
MPI accumulation, zero-weight reduction proof, or long-time validation here.

Next gate: specify a common face reconstruction/flux that accounts for the
residual field, with the approved original-model reduction stated explicitly.
Avoid treating a new second-moment system as already justified by conservation
alone. Preserve existing original-solver trajectories for physical comparisons.
