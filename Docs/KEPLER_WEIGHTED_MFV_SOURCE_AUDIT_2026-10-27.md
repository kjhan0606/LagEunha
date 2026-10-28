# Native mass-flux source audit for the weighted Kepler path

## What is already in the source

`src/Exam/exam.c` contains a compile-time `LAGMFV` face-flux block in the
CPU `getAccVoro2D` path (around lines 6447–6474). It forms a full moving-face
mass flux from the HLLC contact speed relative to the face, adds upwind
momentum transport, and adds an enthalpy flux. A separate `LagMFM` KDK
integrator updates `stress.mass_flux_accum` in two half steps (around
10634–10679). This is a useful reference, but it is not the active CUDA RK4
weight run. The CUDA force kernel returns acceleration and `die` but has no
mass-flux output or RK4 mass state. In the current Kepler launcher
(`SEDOV_PHASE1=1`), `die` is interpreted as total-energy rate; the RK4 host
path subtracts a fixed-mass kinetic-energy rate to recover `ie`.

## Why the existing block cannot be enabled as-is

1. It transports mass across the full face slip. In Kepler shear that slip is
   nonzero for ordinary W=0 Voronoi faces, so enabling it violates the
   approved exact W=0, fixed-individual-mass reduction.
2. In the current state contract pressure is `P=(gamma-1)*ie/V`; `ie` is
   extensive internal energy. The MFV comment calls its added `die` term a
   total enthalpy flux, while the inspected LagMFM KDK half steps add `die`
   directly to `ie`. Separately, the active CUDA RK4 phase-1 path interprets
   `die` as `dE/dt` and converts it with `phase1_ie_stage`, whose current
   `dK=m v·a` assumes fixed mass. These paths do not share one demonstrated
   mass-and-energy contract. Before reuse, derive the kinetic-energy change
   from simultaneous mass and momentum transport and audit the conversion;
   otherwise energy can be counted in the wrong reservoir.
3. The block is not stage-coupled to the CUDA RK4 volume-weight solver and
   does not provide the approved W=0 reduction.

## Required CUDA RK4 integration

The candidate should accumulate per-face antisymmetric rates for mass,
momentum and **total** energy at each RK4 stage. Its additional mass transfer
must arise from the weight-induced face-speed difference, so it vanishes
identically when all weights are zero. Update extensive conserved quantities
(`m`, `m v`, total energy) through RK4, then recover velocity and internal
energy with one explicit kinetic-energy conversion. Keep the existing W=0
RK4 path bitwise unchanged. Any spin extension needs a separate angular
momentum and rotational-energy ledger; the current `GFS_LAGUERRE_ROTATION`
option only corrects face-normal velocity.

The first useful gate is a single native GPU RK4 step from the accepted
jittered `t=6.5` checkpoint, comparing W0 reduction and the weighted flux
candidate with pairwise mass/momentum/energy balances, positive state,
entropy sign, density/orbit response and wall time. The fresh weighted
fixed-mass branch now passes its own `t=10` gates, but that result is not a
local moving-face flux test and does not remove the need for this integration.
