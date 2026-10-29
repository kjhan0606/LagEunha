# Thermal64 RK4 final anchor: isolated candidate

In native v2, the thermal64 K4-to-final reset assigns the authoritative double anchor to U64, then casts it into the float `ie` mirror. The final RK4 combination subsequently uses that float mirror as its base. Thus, after a step has created double precision information below the float mirror's resolution, the next final update discards it. The existing Q1 wall runs have thermal64 disabled, so this finding does not retroactively explain those results.

The attached one-line candidate changes only the thermal64 final base to `GFS_T64_ANCHOR`; the direct-U and legacy branches retain their original expressions. An isolated native executable was built on Ada, without changing running binaries. A separately compiled regression extracts the candidate's actual expression and checks thermal64, legacy and direct-U arithmetic with a double anchor distinct from its float mirror. Its receipt and the build hashes are included here.

This is not yet a native evolution test, a convergence result, or a remedy established for the Kepler failure. Bounded wall events still explicitly reject thermal64; their rollback, authoritative-energy checks and wall energy ledger must be integrated and validated before the candidate can be used in that setting. The current H200 timestep comparison and Ada orbit-cap run use the untouched original executable.
