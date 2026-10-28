# Two-cell inviscid spin and conservation ledger prototype

The [C++/CUDA-compatible prototype](gfs_spin_face_exchange_20261027.h)
performs one isolated face exchange between two cells. It is not attached
to the native hydro solver. The caller supplies the actual mass, momentum
and total-energy transfers, face centroid, and the new cell inertias; the
prototype does not invent a friction or viscosity at `Re=∞`.

The state position is the **mass centroid**, not the Voronoi/Laguerre
generator. Spin `S` and moment of inertia `I` are about that centroid. In
the assumed rigid-spin closure the residual internal energy is
`E - |p|²/(2m) - S²/(2I)`; unresolved internal kinetic structure is not
separately represented. The actual solver must compute consistent centroid
and inertia from its cell state and geometry before this expression can be
used as thermal energy.

For mass and momentum transfer `dm, dp` from `i` to `j` through face
centroid `c`, the prototype transfers the first mass moment by `dm*c` and
the angular momentum by `c×dp`. It updates both centroids and the intrinsic
spins so that the pair sum of `C×p+S` stays constant. Equal and opposite
mass, momentum and total-energy changes preserve those pair sums. Supplying
`dm=0` preserves each mass exactly while an off-center momentum transfer
can change spin. Both cells are updated only if mass, inertia and residual
internal energy remain positive and all values are finite.

The [host test](test_gfs_spin_face_exchange_20261027.cpp), compiled with
warnings as errors, passes off-center pressure-like transfer, mass transfer
with changed inertias, first mass moment, linear and angular momentum,
total energy, translation invariance, fixed-mass reduction, and atomic
rejection of invalid mass, energy and nonfinite input. A small CUDA kernel
that calls the same header also compiled with CUDA 13.0.2 for `sm_80`;
the [device execution test](test_gfs_spin_face_exchange_cuda_20261027.cu)
completed `0:0` as A10 Slurm **437453**. Its
[result](KEPLER_SPIN_FACE_LEDGER_GPU_2026-10-27/result.log) passes pair
invariants and failed-update rollback on an actual GPU. The one-thread test
is not a native fluid simulation or a parallel performance measurement.

This is an accounting building block, not a complete spin fluid model.
Full evolution needs a physically specified inviscid face flux and its
energy, a consistent mass-centroid/first-moment volume term, polygon inertia
and shape evolution, generator motion and RK4 staging, and race-free GPU/MPI
face accumulation. The existing variational generator force alone does not
uniquely determine the spin source, as shown by the
[pair-torque analysis](KEPLER_VARIATIONAL_SPIN_PAIR_TORQUE_2026-10-27/README.md).
The approved W0 fixed-mass reduction and local/long-time Kepler tests remain
mandatory before any native promotion.
