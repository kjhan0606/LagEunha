# Weight-motion mass exchange prototype

The isolated [weight-motion helper](gfs_weight_motion_remap_20261027.h)
constructs a pair transfer to the existing spin/affine energy ledger. Let
`n` point from site i to j, `d` be the site distance, `ΔW=W_i-W_j`, and
`Δv_n=n·(v_j-v_i)`. Relative to the Voronoi face at the same site positions,
the exact instantaneous power-plane speed correction is
`ΔWdot/(2d)-ΔW Δv_n/(2d²)`. The first term is weight change; the second is the
finite-weight correction as site separation changes. Positive speed expands
i and transfers mass from j to i. The helper uses the swept donor's density,
velocity, specific total hydro energy, and nonrigid energy; it does not
derive those quantities. Momentum and energy are transported consistently
with the donor mass and each pair transfer is equal and opposite.

When all weights remain zero, their rates are zero and this weight-only
transfer vanishes exactly, consistent with the approved fixed-mass W0
reduction. A zero weight at one instant is insufficient if its rate is not
zero. Equal weight rates also have no geometric effect.

The host test includes a uniform moving-face case, weight-rate and site-
separation contributions, the W=0 reduction, both donor directions, failure
rollback, and 50,000 deterministic face reversals and Galilean boosts. It
checks algebraic symmetries, not the accuracy of an evolved disk.
On 2026-10-27 it passed with GCC using `-std=c++14 -O2 -Wall -Wextra -Werror`.
CPU-only Slurm job **437563** verified the earlier staged input hashes and
passed CUDA translation and link. Its GPU job **437562** passed the original
weight-rate-only test. After adding the site-separation term, A10 job
**437598** failed a test expectation because that case accidentally used
`dt=1` while expecting `dt=0.1`; it did not show a kernel or ledger defect.
The corrected CUDA test **437600** passed on an A10. A corrected test for an
actual evolving hydro step remains future work.

This is an isolated geometric remap of frozen face states. Its upwind
direction follows the weight-induced face-speed correction; that is generally
not the upwind direction of the full ALE/Riemann mass flux. Combining split substeps can
therefore introduce errors. Pressure traction and work, gravity, evolution
of polygon mass moments and inertia, RK4 stages, MPI and restart integration
remain undefined. It is not integrated into the user's production solver.

The earlier [native one-step audit](KEPLER_WEIGHT_ONLY_MASS_REMAP_2026-10-27/README.md)
found that a first-order weight-only mass remap almost cancels the apparent
fixed-mass density improvement. This finding must be revisited in an actual
conservative evolution; neither that one-step estimate nor this prototype
proves long-time accuracy. The [W0 consistency analysis](KEPLER_W0_ALE_FIXED_MASS_CONTRACT_2026-10-27.md)
explains why exact moving-face ALE flux and fixed Voronoi-cell mass are not
generally the same constraint.
