# Same-IC fresh weighted Kepler at t=6.5

The 5%-jitter, bounded 55,351-site native calculation **437499** reached
exactly `t=6.5` from the accepted `t=0.825` restart, using v5's CUDA+MPI+
OpenMP RK4 solver, one MPI rank, no steering, fixed cell masses and evolving
weights. Its [weighted receipt](weighted_receipt.json) passes all 21 gates,
including equilibrium density below `1e-3`, independent power areas, zero
energy injection, energy, orbit and positivity. The matching original W0
checkpoint passes its [19 checks](w0_receipt.json).

The [independent comparison](comparison.json) confirms the same initial hash
and exact endpoint time. Disk mass-weighted mean absolute density change is
`1.69828e-3` W0 versus `5.01080e-5` weighted, mean absolute log-entropy change
is `1.49389e-5` versus `6.91683e-6`, and relative total-energy errors are
`5.24739e-9` versus `5.24667e-9`. The [density-field figure](density.png)
uses the actual Voronoi/power winner at each pixel and identical colors. The
lower panels show each cell's density relative to that same cell's initial
density; their means are mass weighted. The two runs have different binaries
and restart trajectories. This is a fixed-mass experimental weighted method;
cell density responds to weighted geometry and local moving-face mass/energy
consistency has not been demonstrated.

CPU job **437545** postprocessed saved velocity gradients and uniform-density
polygon moments at this endpoint: [W0](w0.json) and [weighted](weighted.json)
rigid-rotation projections are about 11.33% of reconstructed affine kinetic
energy. This is not a native evolved spin and not irreversible heating. If
the projection were subtracted from existing thermal energy while artificially
holding the old total-energy reference fixed, approximately 15.6% of disk
mass would have nonpositive residual thermal energy; that counterfactual is
not an initialization prescription. A physically consistent spin model
requires a compatible energy reference and torque/work closure.

The gated continuation **437530** starts from this exact accepted checkpoint
and runs to `t=10`, with **437531** conditional on its success to `t=20`.
The `t=6.5` density result is not the full 384² periodic control, proof of
physical local face transport, or the GIZMO `t=600` benchmark.
