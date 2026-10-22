# Resident mass and full positive-rule construction — 2026-10-22

The original long Kepler goal is ACTIVE. The accepted trajectory remains at
**t=8.2355880701376864e-5, four accepted steps**. Neither t6–10 artificial
heating suppression nor the longest GIZMO comparison is established.

## Actual next high-rule RHS: adopted and retained

`material_full_resident_rhs_v1_ada_mpi4` finished native0 and runner0.
Original native MPI entry, MPI4/OpenMP8, ONE RTX5000 Ada; unchanged71080cells,
1138432triangles,10816441scalar DOFs, and2040256512positive points/pass.

New inner PCG x/r/z/p/Ap remain DD on GPU, including owned-only DD dot and AXPY.
Sparse ghost hi/lo pairs use pinned staging and canonical peer-order scatter.
There are no unordered double atomics or whole-vector transfers per iteration.
The authoritative outer residual is still computed with the original saved
binary80 M; zero-RHS handling,256inner and3outer bounds, and1e-16 tolerance
remain. This does not claim DD and binary80 are identical arithmetic.

| Actual same-RHS measurement, Ada | Value |
|---|---:|
| New resident two-component mass solve | 16.3184886s |
| Previous host-vector DD reference | 30.0907918s |
| True residual x/y | 8.12765e-19 / 8.12574e-19 |
| Inner iterations x/y; outer corrections | 54/54;1/1 |
| Reference work difference/U | 3.99881e-19 |
| Reference torque difference/(2K) | 3.80824e-23 |
| Worst full-cell E rate error/U | 3.80327e-13 |
| Worst full-cell L rate error/(2K) | 4.01116e-16 |

The recorded halo wall includes device synchronization/product completion;
14.67s is not a pure MPI communication measurement. The actual full-vector DD
traffic was173146560bytes/rank across both component corrections, plus sparse
ghost traffic27070848bytes/rank. The comparison was performed once on the actual
next RHS, not by replaying an accepted step.

`material_full_resident_rhs_audit_v1/contract.json` independently checks74hashes,
all71080physical/rotation ledgers, and every saved k1 xdot against the accepted
velocity **bits**. It passes. Prepared next high-rule k1 is retained at
`material_full_resident_rhs_v1_ada_mpi4/run/next_RK_stage1/prepared_k1_r[0-3].bin`.
The payload binds accepted coefficient/immutable manifest hashes, time, rate
ledgers and actual acceleration. Reuse it for the same high-rule first stage.
Do not recompute its force/mass solve merely to initialize a new executable.

## Positive support feasibility: complete representative evidence

`material_positive_support_v1` used one fixed nonnegative least-squares support
proposal on all24retained matrices, with every324–400constraint row present.
No singular direction, density tail or physical state was deleted. The package
signature/behavior was checked against installed SciPy1.18 and its
[official nnls documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html).
This Python operation proposes integration weights; physical fields and
evolution remain native C++/CUDA.

`material_positive_rule_audit_v1_mpi4` completed native0 in1.67179213sec:
**24/24 pass**,93184original points ->5076strictly positive original-subset points
(~18.36x for this deliberately mixed representative set). All-row residuals are
evaluated independently in binary80, then contracted with the actual accepted
x/v, saved actual next acceleration and unchanged saved M. Tests include cold-U
normalized K/Kdot/heat, angular momentum, corrected-J/dJdt/gradient, triangle-mean
pressure/volume-force proxy and a finite arithmetic QJ bound. Centered x/v
match `st_point`'s origin subtraction; using absolute coordinates would give an
irrelevant cancellation bound. Every candidate retains original generator/source
grouping, fixed M and generating entropy. Full-cell force and full RK-stage
oracle comparisons remain required. Numerical table rank is not certified or
used to silently remove constraints; candidate correctness is checked on ALL
rows and actual-state contractions.

This is a positive finite cubature construction rather than a Gauss-order
sweep. No reduced production operator or time step has been committed.

## Full-domain native construction: live

`material_full_positive_rule_v1_ada_mpi4` builds/runs native C++/CUDA/MPI4 from
the retained accepted state and prepared k1. Its four GPU ranks stream the same
original reference points, source-complete geometry/thermal contexts and DD
constraint rows. Sixteen private local positive-support workers use four BLAS
threads each. Native C++ validates every proposed rule in binary80, preserving
the original positive rule on candidate failure. The whole-domain geometry pass
is needed for source contexts; the24saved representative matrices/supports are
reused, not rebuilt/reoptimized as separate tests.

No full constraint atlas is written:32MiB GPU matrix output/rank,1048576point
cache/rank, four bounded jobs/rank. Only compact rules, signed moment residuals,
rounding bounds and discrepancy/fallback ledgers are retained. Planned compact
output <=24GiB and worst bounded host budget160GiB. Available host memory was
954GiB and free GPFS space216TiB at launch; GPU was idle. Current field/time/M
and entropy stay unchanged throughout preprocessing.

Runner tool **9211** is live as last verified. Supervisor/native PID files and
`worker_pids.json` are in this directory. Actual native MPI launcher PID2166274;
GPU ranks2166305–2166308. Do not duplicate or restart on an observation timeout.
At elapsed~14m46s, atomic committed prefixes contain12672/12096/12644/13204
triangles (50616total,~4.4%), all compressed:90721792points->10996589 (~8.25x).
Projected total ~5–6h. Native launcher/fourrank PIDs were verified live. This is
an early cost estimate, not completion. The producer stops with retained prefix
if projected preprocessing cost exceeds24h after2048triangles/rank, or at its
bounded wall budget. Candidate-ready also requires every triangle, original
2040256512points, total reduction>=4x and cold-U normalized anchor U/Phi bounds.

## Remaining connection before long production

Read the actual terminal result; an exit0 alone can mean a safe partial/cost
stop. Validate full compact rule/residual records and global discrepancy budget.
Connect these rules to evolved/cache native point selection, with full-rule
fallbacks, immutable generation and authoritative M unchanged. Compare the
actual next RK step at the same stage states against the original rule, reusing
the saved high k1. Preserve additional quadrature error/work ledgers and original
dt bounds. Then assess measured step cost before the t1/24h trajectory.

Original rule and all accepted data stay authoritative. Sampled positivity is
not global map injectivity; the enriched-space spectral CFL, material topology
lifetime, actual conservative repartition/recipient thermal compatibility,
W0 dynamic Voronoi reduction, shocks and paper-matched boundary/normalization
are unfinished. Moving constant-Pext material BC differs from original periodic.
No heat correction, damping, steering or separate spin reservoir was added.

## Native connection continuation

See RESULT_POSITIVE_RULE_CONNECTION_20261022.md for the compiled paired-step
connection and live single dependency PID2296727. The full producer is unchanged
and running. No new accepted step is claimed until the actual downstream result.
