# Midpoint-offset pressure eligibility

See ../../Docs/UNIFIED_HYDRO_MIDPOINT_PRESSURE_2026-11-05.md for definitions,
results, independent 1D derivation, and limitations. The current common face
fit and original native geometry define a midpoint-offset face divergence D.
The exact transpose includes the gradient-fit transpose. Algebra passes, but
cell-average linear-pressure force-density errors do not decrease under a
matched factor-two refinement. Nodal pressure also fails in 2D/3D. The candidate
is rejected at the preregistered pressure gate before Bloch/native stiffness,
weight solves, or any PDE evolution.

Original native small geometries are reused and exactly tiled. One complete
central primitive pattern is tested with a two-layer seam guard. Archived
native doubled geometry verifies the tiling scale. No new multidimensional
mesh, global coefficient solve, or fluid step is performed.

The 1D W=0 D=B reduction is preserved. The inherited cell-average pressure
error is separately derived and verified with a C interval proof and the
existing native variational-pressure helper. The earlier constant-pressure
contact initialization cause remains valid. A point quadrature model and a
literal cell-mean model must not be silently interchanged.

historical_reference contains the Nov2 formula for comparison, not a runnable
reproduction or selected solver. Its original fit differs from the common
fit used in the new matched-formula comparison. Source snapshots, commands,
outputs and SHA256 provenance are retained. Original production state and
manuscript benchmark remain unchanged. The full objective remains incomplete.
