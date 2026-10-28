# Existing general HLLC helper parity and active-path distinction

Use10736 cell-constant donor-state samples at certified segment midpoints plus
10 manufactured controls. Actual samples use phase1=0 for comparison; controls
exercise phase1=0/1, extreme pressure ratios, expansion/compression and contact.
Gamma5/3, postype double. These are not the fitted density face reconstruction
or a replay of all active production caller options.

Extracted existing CPU/CUDA helpers from the prior provenance-locked probe,
preserving arithmetic. A thread-local phase1 selector enables safe per-input
CPU controls. omp_set_dynamic(0), actual team sizes measured as1 and4.
Same GNU13/CUDA13 FMA-off compilation; actual CUDA kernel executed.

-10746 cases finite; CPU1/OpenMP4/CUDA bitwise identical, max difference0.
-Uniform boost7: star-P error<=6.47e-17, star-velocity shift error<=1.78e-15.
-Reverse normal and side ordering: P error0, velocity error<=1.12e-16.
-Pressure-ratio control P_L=1,P_R=1000 returns P*=1 in phase1=0 guard,
  versus499.512 in phase1=1 HLL-average branch. This confirms different
  branch behavior, not correctness against an exact general Riemann solution.

## Do not misattribute the generic expansion guard

The generic hllc_face_2d branch tests uR-uL>max(cL,cR) and returns tiny pressure.
788 actual comparison samples satisfy that predicate;308 are below the
ideal-gas vacuum threshold2(cL+cR)/(gamma-1) used explicitly by the separate
existing tcf_pstar routine. These counts are sample evaluations, not distinct
physical shocks, nor proof that the generic guard alone is a production bug.
No Roe/last-resort fallback coverage is claimed from these inputs.

Read actual running process environments for PIDs2286374 and2295268, limited
to relevant solver keys. Both have GFS_TCF=1,GFS_TCF_QBP=1,GFS_TCF_FACE_Z=1,
GFS_QBP_CHI=1,SEDOV_PHASE1=1; HG flags are also active. Their params.dat uses
av_mode5. Therefore the generic helper probe is not a complete description of
the active force/pressure calculation. The distinct TCF path must be inspected
and tested before attributing its effects to the generic HLLC guard.

## Next priority

Trace actual TCF caller selection and compare its CPU/CUDA pressure evaluation
on the same inputs, including the intermediate expansion states. Retain the
old generic behavior and production binaries until evidence identifies a
relevant defect. Backend agreement alone does not prove physical accuracy.
No source correction, restart, MPI validation or long-time stability achieved.

Evidence: `canonical_existing_interface_backend_2026-10-18.json`,
`canonical_existing_interface_transforms_2026-10-18.json`.
Scripts: `prepare_existing_interface_probe.py`, `build_existing_interface_probe.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_existing_interface_backend_20261018`.
