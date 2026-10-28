# Active TCF is a CPU force path; generic CUDA parity is not TCF parity

Authoritative current source and live configuration checks:

- `Exam/exam_gpu_extract.c`, getAccVoro2DBlend_GPU explicitly aborts on
  GFS_TCF or GFS_DYN_W2, explaining GPU force equivalents are not implemented.
- No tcf_pstar/TCF force branch was found in the scoped exam_gpu.cu search.
- Active fixed-HG params.dat has GAS gpu_enabled=0. Both inspected baseline
  and fixed-HG process environments have OMP_NUM_THREADS=8 and GFS_TCF=1.
- CPU caller requires one face point, SEDOV_PHASE1 total energy, no schemeH
  or dynamic weights, and equal face weights within1e-12*d^2. It therefore
  is not already a general moving-weight Laguerre solver.

Consequently earlier generic HLLC CPU/CUDA helper parity does NOT validate
the active TCF force path on CUDA. Do not describe these active TCF trajectories
as GPU-force runs. This does not imply all CUDA geometry/other modes were
removed, and no original implementation was deleted or replaced this turn.

## Actual CPU TCF helper comparison

Extracted the existing tcf_pstar arithmetic unchanged into a small CPU library,
recording source hash. Evaluated the same10736 raw cell-constant donor samples:
9940 positive outputs and796 zero outputs. Of308 samples where the generic
large-expansion guard triggers below the ideal-gas vacuum threshold, only10
have positive TCF pressure. Thus the earlier308 count is not308 demonstrated
pressure disagreements. Exactly10 actual outputs differ by more than1e-12;
maximum absolute difference4.42501e-6. Others agree to the tested threshold.

This is not an exact Riemann accuracy test or a replay of reconstructed active
faces. Furthermore the TCF caller uses the pressure function in a difference
quotient to form effective impedance Z, then downstream QBP/force logic; the
returned P* is not by itself the full applied force. No causal claim about
the current Kepler failure follows from these10 differences.

## Priority correction

Separate two required workstreams explicitly:
1. Stabilize and validate the actual CPU/OpenMP TCF physics against long Kepler
   trajectories, preserving the historical t6--10 and t20+ acceptance scope.
2. Implement and validate the SAME full equations on CUDA and MPI before
   claiming the user's cross-backend objective. Porting only tcf_pstar or passing
   generic HLLC tests is insufficient; geometry, Z/QBP, work and reductions matter.

Do not enable GPU on TCF merely to bypass the unsupported guard. Do not treat
the experimental remap studies as an existing Laguerre extension of TCF.
Next trace the actual impedance/QBP assembly and identify its complete data
dependencies and existing tests before deciding a scoped port or physics change.

Initial analysis command used an incorrect relative sibling path and failed
before producing results; corrected analyzer uses an absolute root. The small
library was not rebuilt or altered to obtain the final comparison.

Evidence: `canonical_active_tcf_pressure_2026-10-18.json`, source/hash snapshot
`tcf_pressure_audit_2026-10-18.c`, `tcf_pressure_provenance_2026-10-18.json`.
Analyzer: `analyze_active_tcf_pressure.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_active_tcf_pressure_20261018`.
