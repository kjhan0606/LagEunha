# Bounded Q1 zero-time thermal64 initialization at t=10

The existing native v2 executable initialized its thermal64 payload from the preserved, accepted t=10 bounded Q1 checkpoint on the Ada host. `GFS_THERMAL64_INIT_ONLY=1` stopped the process immediately after initialization: no physical time step and no wall-event evolution occurred. The output is a diagnostic snapshot, not an accepted restart. The native wall-event/thermal64 compatibility guard was not modified or bypassed for evolution.

An independent all-site check passed for 90,352 noncontiguous global IDs. Position, velocity, mass, weight and the original float internal energy were unchanged; every double U was the exact promotion of that energy, the tags and anchors matched, and float mirrors and EOS checks passed. Independent clipped Voronoi reconstructions with two guard radii agreed with each other to 1.76e-10 maximum relative error. The promoted double volume agreed with the independent geometry to 1.20e-10 maximum relative error; the largest logarithmic volume change relative to the previous float volume was 5.95e-8. This geometric precision change is not physical heat.

Source checkpoint SHA-256: `58afbe3ab219ddef0b8dbb9b032789236ec01815bcf80457a5d73e6b90bd6c28`. Diagnostic snapshot SHA-256: `0456a8c6dafb3f7ee7e681089c8c475462a4fe7e4ee062a19cf19f4f2db52381`. The complete native log and snapshot remain in `/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/thermal64_init_t10_v1/`; the audit receipt is included here.

This experiment neither repairs the inherited core error nor demonstrates stable or conservative thermal64 evolution with bounded wall events. That combination requires a separate implementation and validation of the full rollback, wall work and energy ledger.
