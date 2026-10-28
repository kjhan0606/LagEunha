# Two-pass polygon centroid and inertia for a future spin state

The [C++/CUDA-compatible routine](../gfs_polygon_moment_20261027.h) obtains
area, centroid offset and central polar moment from consistently oriented
polygon edges. It makes two linear passes over the faces: the first computes
area and centroid by boundary sums, and the second computes the moment in
centroid-relative coordinates to reduce cancellation. With uniform cell
density, the scalar 2D inertia is `I=m J_c/A`. The generator need not be at
the centroid. No subcell quadrature or steering is used.

The [host test](../test_gfs_polygon_moment_20261027.cpp) passes analytical
square, offset rectangle, triangle, reversed orientation, distant-generator
and invalid-polygon cases. An [independent comparison](../audit_kepler_polygon_moment_native_faces_20261027.py)
then rebuilt faces from hash-verified native snapshots and compared the
same routine against triangle-by-triangle area, centroid and inertia
integration:

| Checkpoint | Disk cells | Directed disk faces | Maximum inertia relative difference | Maximum area difference from native snapshot |
|---|---:|---:|---:|---:|
| [W0 jitter t=6.5](w0_jitter_t65.json) | 49,038 | 294,228 | `8.88e-16` | `5.94e-8` |
| [Weighted jitter t=10](weighted_jitter_t10.json) | 49,038 | 294,228 | `8.88e-16` | `5.96e-8` |

The largest centroid difference against triangle integration is below
`1.74e-18` in these coordinate units. Both calculations use independently
reconstructed faces. The weighted `t=10` snapshot failed the
equilibrium-density gate and is used here only as a geometry diagnostic.

The [face-capture audit](../audit_kepler_polygon_moment_gpu_face_capture_20261027.py)
also uses actual CUDA face endpoints saved during the earlier native W0
`t=6.5` step. For 64 [aligned](gpu_capture_aligned_t65.json) and 64
[jitter](gpu_capture_jitter_t65.json) cells, the routine agrees with triangle
integration in inertia to `4.44e-16` relative; its areas agree with the
native force-input volumes within `2.95e-8` and `4.80e-8`, respectively.
The moment routine itself ran on the CPU in this capture check; the full
weighted native GPU face arrays have not yet been checked this way.

The [synthetic CUDA batch test](../test_gfs_polygon_moment_cuda_20261027.cu)
completed `0:0` as A10 job **437458**. For 55,351 four-face rectangles over
100 timed launches, the [result](gpu_synthetic_result.log) reports mean
kernel time `0.0757 ms`, no invalid cells and maximum inertia relative error
`5.56e-16` against the analytical rectangles. This synthetic kernel timing
excludes native tessellation, host/device transfers, MPI and real polygon
face-count variation; it is not an end-to-end speed measurement. Native
FaceCSR attachment, realistic GPU cost, nonuniform subcell density, 3D
inertia tensors, spin flux and the long-time Kepler test remain open.
