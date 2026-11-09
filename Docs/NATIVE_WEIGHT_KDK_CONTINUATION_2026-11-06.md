# Native weighted KDK continuation and retry reduction

Artifacts: Research/MaterialVolume_20261105/native_weight_bins_20261106_v1.
All runs use the same A10 native solver, one process, two OpenMP threads.
No original IC replay was used for this continuation: each case starts from
its saved t=.02 native-dt checkpoint. Per-segment lineage includes source path
and SHA256. Accepted endpoints are .035,.05,.065,.08,.095,.1.

## Time selection

Within a continuous process, the first computed native CFL is the bin root.
Select the largest root/2^k not exceeding the current bound, then apply the
endpoint ceiling. Later-stage bound and full-payload retry remain unchanged.
No CFL coefficient or physical parameter was increased. This is global
dyadic stepping, not individual time stepping. The root resets on restart;
exact uninterrupted/restarted trajectory identity is not asserted.

Control and binned runs share the same t=.02 checkpoint and t=.035 endpoint:
- Gresho: retries 28 -> 1, steps 29 -> 30; summed step wall time
  1.525101 -> .953726 seconds.
- KH: retries 37 -> 1, steps 38 -> 39; summed step wall time
  1.641680 -> 1.039864 seconds.
This approximately 37% reduction is one small internal experiment, including
initial GPU work in step timing. It is neither a repeated throughput benchmark
nor a GIZMO runtime comparison.

## t=.1 result

Both checkpoint chains reach .1 with exactly unchanged per-cell mass.
Relative augmented total-energy errors against their original IC:
Gresho 6.127e-10; KH -7.268e-11.
Maximum adiabatic entropy change is below 4e-15.
The absence of entropy production is intentional for this smooth branch and
does not qualify any shock test.

Weight kinetic energy at .1 is .00129514 in Gresho, approximately 1.59% of
translation kinetic energy; KH weight kinetic energy is 1.93308e-6,
approximately 1.43e-6 times translation kinetic energy. Conservation alone
does not justify that auxiliary energy reservoir physically. That question
and comparison with the original W0 method remain required.

## Incomplete requirements

The grids are Gresho 16x16 and KH 16x32; these are continued low-resolution
engineering tests, not a resolved instability result. KH/RT physical
endpoints, shock/entropy consistency with the successful 1D material Riemann
solver, general MPI, and the <=2x GIZMO end-to-end target remain unverified.
3D is deferred until 2D completion, as instructed.
