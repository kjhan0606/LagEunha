# Actual periodic 2D evolution with the common operator

See `../../Docs/UNIFIED_HYDRO_PERIODIC_2D_2026-11-04.md` for methods, measured
results, the active-set failure and correction, and unresolved limitations.
All hydrodynamic evolution uses the original C/OpenMP+MPI+CUDA executable.
Python scripts orchestrate runs and analyse saved outputs.

The first contact/shock batch uses general2d_build_v1. Two shear runs failed
with a nonzero exit code. The corrected shear batch uses general2d_build_v3
and includes a third resolution. Both failure and success evidence is retained.
The captured graph replay isolates sub-ulp active-set cycling. The correction
adds no constitutive switch or tolerance. The 1D nine-case regression produces
byte-identical final fields after the correction.

Native binary snapshots, static archives, objects, executables, and the
captured failing graph remain under the absolute GPFS paths in manifest.json.
The manifest records hashes of the retained large data and publication files.
Build command records are site-specific hybrid relinks, not portable clean
rebuild instructions. The current source snapshot is distinguished from the
earlier v1 binary. Run scripts refuse to overwrite existing run directories.

This archive does not establish complete shock qualification, 3D evolution,
MPI scaling, long Kepler stability, or a matched GIZMO performance result.
The axis-aligned Sod transverse velocity requires a separate causal audit.
The manuscript remains the previously published baseline snapshot until the
experimental common scheme is qualified.
