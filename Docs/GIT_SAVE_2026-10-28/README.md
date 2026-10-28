# Source and research preservation, 2026-10-28

The legacy September checkout is preserved on archive/lagforce-src-20261028, and the October 17 validation work on codex/physics-validation-20261017. Their exact commit and push receipts are alongside this file.

This research branch preserves the actual 42-file native baseline source overlay, the accepted clip-swap source, restart and physics audit records, density figures, the verified Slurm cleanup and the latest default-off GFS_REUSE_TESS_PADDING experiment. The latter changes only redundant host halo reconstruction while retaining post-gradient synchronization and index checks; its native comparison job is 439001. It is not enabled in the running long-time trajectory.

The two-A100 MPI gate at t=1.55 passed with a measured 1.08714345x speedup; this does not establish long-time accuracy. The native long-time and refined-resolution jobs continue. GIZMO t=600 remains unverified. Native binary snapshots and compiler outputs remain on GPFS.
