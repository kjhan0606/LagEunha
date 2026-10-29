# Frozen Q2 affine weight exchange check

The actual `gfs_weight_exchange_segment_affine` C candidate is tested on the
bounded Q2 t=47.50031653561342 geometry. The 3,931 inner-disk cells have
12,508 touching faces; synthetic signed power weights and weight rates are
used with fixed generator positions. A shared affine velocity field, unit
uniform density and constant specific internal energy specify the transport.
An independent power-diagram reconstruction computes polygon integrals at
two perturbed weights, whose centered difference is compared with the C face
operator. The [audit](KEPLER_WEIGHT_AFFINE_GEOMETRY_2026-10-29.json) checks
all 13 raw mass, spatial, momentum, momentum-position and total-energy
moments, as well as equal-and-opposite shared-face contributions.

For these assumptions, weight-driven mass transport changes the geometric
moments by exactly the same face integral as the weight-driven flux. Thus the
weight contribution alone cannot fix the pre-existing moment-rate mismatch
between fixed-mass Voronoi geometry and material affine motion: at this
uniform-density frozen-state limit it cancels on both sides of that comparison.
This is not a general impossibility result for nonuniform reconstructed
states, coupled weight laws or an evolved variable-mass model.

This test establishes face integration and signs for a chosen smooth field,
not physical closure, moving-generator flux, pressure work, positivity,
native negative-weight support, MPI, or finite-time spin/Kepler stability.
The candidate remains disconnected from the native simulation.
