# Dimensional conservative FV primitives: C/C++ and actual Ada checks

This records the primitive-only checkpoint. The subsequent native CPU
evolution implementation and its limits are reported separately in
`UNIFIED_HYDRO_FV_NATIVE_2026-11-05.md`.

This implements the first native pieces of the grouped Astra plan in
`UNIFIED_HYDRO_FV_ALE_PLAN_2026-11-05.md`. It does not release a new fluid
solver. `Exam/gfs_fv_ale.h` is included only by its native audit programs;
existing production pressure, geometry, viscosity, integration and restart
layouts remain unchanged.

## One shared mathematical definition

The header provides primitive recovery from the declared FV mean-state
energy, one frozen lab-frame HLLC/HLL fan, complete sampling `F-w*Q`, and
the first admissible forward-Euler ray interval. The latter finds the first
zero of m or `2mE-|p|^2`, using the stable quadratic roots. It returns a
negative status for invalid input and infinity for no finite boundary.
Callers must still apply the registered safety factor, validate actual
stage states and actual geometric V, and enforce CFL/GCL checks. No floor,
thermal reservoir, kinetic subtraction, or endpoint state repair is used.

The final storage is always `[m,px,py,pz,E]`, with unused momentum slots
EXACTLY zero. There are still only d+2 independent variables. The fixed
layout makes energy location and payload component bounds identical in C
and CUDA and avoids the failed runtime packed-component path below. This
is an experimental internal payload, not a change to accepted restart ABI.

HLLC uses the existing Davis primary bounds and positive star-state check.
An inadmissible HLLC state selects the COMPLETE HLL flux with symmetric
local-LF outer bounds in the normal mean frame. The HLL intermediate state
is checked; no residual or mass flux is removed.

## Independent native C/C++ evidence

`Exam/HydroCore/fv_flux_audit.c` checks 40 deterministic fan fixtures in
each dimension, including a strong expansion forcing HLL fallback. It
checks face reversal, Galilean transformation of the COMPLETE energy and
momentum flux, independent uniform-state physical ALE flux, and the mass
derivative `dFm/dwn=-rho` inside each sampled fan region.

The same planar Riemann data embedded in 1D/2D/3D produce identical active
components and exactly zero inactive momentum flux. Five independently
specified cone boundaries per dimension include linear/quadratic roots,
a double boundary, both positive roots, and a boost. States on the fixed
0.9 interior fraction are positive, states just beyond the first boundary
are rejected. Invalid/zero-rate controls are separate. This is not a
claim that the future PDE timestep has already been qualified.

In the C++ build, all nonfallback 2D flux components also match the existing
`gfs_ale_fan.h` implementation, with energy mapped to the fixed final slot.

| Check | Final maximum normalized discrepancy |
|---|---:|
| Face reversal | 6.15e-16 |
| Galilean complete flux | 3.86e-15 |
| Uniform physical ALE flux | 3.54e-16 |
| Finite-difference mass derivative | 6.71e-11 |
| Dimensional embedding | 0 |
| Analytic cone roots / boost | 1.61e-15 |
| Original 2D nonfallback HLLC | 0 |

Both optimized native C and C++ processes exit zero with the preregistered
limits. These fixtures exercise the new shared functions rather than
replacing actual native geometry or evolution with Python.

## A failed CUDA implementation was diagnosed before any coupling

The initial packed layout put E at `dim+1` and traversed `dim+2` components
in local fan aggregates. Its optimized CUDA audit failed: cone intervals
and intermediate state diagnostics differed substantially from CPU. The
same source compiled with device debug optimization disabled agreed to
9.29e-17. Compute Sanitizer reported zero memory errors on the failing
binary. This does not by itself establish the underlying compiler defect.

Small isolated call combinations preserved the original state correctly
when dimension was a compile-time constant, but runtime dimension plus fan
sampling could change subsequent scalar diagnostics even when copied raw
state values agreed. Forcing inlining did not fix it. Constant component
loads in thermal recovery alone did not fix the complete audit either.
All failed binaries, source snapshots and controls are retained. No PDE,
accepted checkpoint, or physical model was adjusted to hide these failures.

The final fixed five-slot payload and fixed component traversal remove the
runtime packed energy remapping. The identical physical C/C++ checks still
pass, and the optimized CUDA audit now passes all 96 fixtures, including
28 HLL fallback cases, on the actual RTX 5000 Ada:

    Maximum CPU/GPU normalized discrepancy = 3.2058011602380016e-15.
    Compute Sanitizer: 0 errors; application exit 0.

The CUDA audit retains mixed runtime dimensions, rather than avoiding the
failed path by testing only a single known dimension. Raw CPU/GPU traces
include the fan flux, pressure/cone diagnostics and normalization data.
The passing build uses ordinary `-O3`, not debug mode, relaxed physics
checks or a larger tolerance. No general claim about other compilers or
every possible fan state is implied by this finite audit.

## Authoritative scope and next implementation

The first fixed-layout pass is v6. The final run, including explicit
rejection of nonzero inactive momentum rates, is
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/fv_primitives_v7`.
Versions v1--v5 preserve the failed implementations and diagnostics;
v3 includes debug-build and memory-check controls and v4 the isolated
call combinations. Compact evidence is in `Research/FVPrimitives_20261105`.

Native geometry adapters, centroid reconstruction and fixed limiter,
screened weight solve, restartable birth mu/m/W, SSPRK2 and rollback/GCL
remain to be implemented and tested under the SAME grouped plan. No FV
PDE step or matched speed measurement exists yet. The separate corrected
Shu trajectories retain the previous W0 variational scheme. Full common
1D/2D/3D physical accuracy, the remaining 2D suite and <=2x GIZMO remain
unfinished, and the complete user goal is active.
