# Exactly uniform pressure in the experimental weight-rate projection

The experimental CUDA weight solver projects its rate `r` onto `h·r=0`,
where `h` is the pressure-weight work gradient. Its v3 host code aborted
when `h·h=0`. If `h=0`, however, the constraint is already satisfied by
every finite rate; the least-change answer is the original rate with
projection coefficient zero. An exactly uniform cell pressure makes every
pressure-difference contribution to `h` zero.

The accepted fresh jitter `t=0` snapshot stores the same float pressure
`9.999999974752427e-7` in all 55,351 cells. Its actual native stage logs
show nonzero projection coefficients, so that run did not demonstrate an
exact-zero runtime failure; stage pressure/geometry can differ from the
stored initial field. The source-level zero-gradient branch is nevertheless
unambiguously reachable for an exactly uniform pressure input.

The isolated [projection helper](gfs_pressure_neutral_projection_20261027.h)
leaves rates unchanged and returns coefficient zero for an active, finite,
exact-zero gradient. It retains the original summation and update order for
nonzero gradients, keeps an empty active mask and nonfinite inputs as errors,
and rejects a nonfinite coefficient. No arbitrary cutoff is applied to a
small nonzero gradient. The [host test](test_gfs_pressure_neutral_projection_20261027.cpp)
passes zero, ordinary nonzero, `1e-30`-scaled, inactive and nonfinite cases;
the ordinary nonzero output matches the previous arithmetic bitwise. A
later local host test additionally confirms that an overflowing coefficient
is rejected without modifying the rate.

The [isolated source patcher](apply_kepler_zero_pressure_projection_20261027.py)
was smoke-tested against the exact v3 source. CPU Slurm **437448** built the
patched CUDA binary `0:0` and repeated the host test. Its binary SHA-256 is
`3120aea6dad944b144e04d4597cb3f644abc6c803ff2318b5862b6bf8eec0adc`.
The [fresh native one-step runner](kepler_q2_zero_pressure_v4_t001_20261027.slurm)
was submitted as **437452** to compare its physical fields with the
completed v3 `t=0→0.001` weighted pilot. It completed `0:0` in ten seconds;
the [candidate audit](KEPLER_ZERO_PRESSURE_PROJECTION_2026-10-27/v4_audit_t001.json)
passes all 18 checks and the [paired comparison](KEPLER_ZERO_PRESSURE_PROJECTION_2026-10-27/v3_v4_pair.json)
finds bitwise equality in all ten initial and endpoint physical arrays.
The input and executable hashes are recorded in the receipts. That native
run exercises the ordinary nonzero-gradient path. The existing v3 long jobs
retain their original binaries. This numerical repair does not establish
local moving-face conservation or long-time Kepler stability.

The separate [two-cell CUDA regression](test_kepler_zero_pressure_gradient_cuda_20261027.cu)
uses the actual native face-CSR Jacobi and pressure-gradient probe with
equal pressures and nonzero raw rates. A10 job **437455** completed `0:0`:
the [old result](KEPLER_ZERO_PRESSURE_PROJECTION_2026-10-27/gpu_old_result.log)
shows zero gradients and status `-3`; the
[fixed result](KEPLER_ZERO_PRESSURE_PROJECTION_2026-10-27/gpu_fixed_result.log)
shows status `0`, coefficient `0`, and bitwise unchanged nonzero raw rates
`(+0.1999999936,-0.1999999936)`. Both use the actual CUDA face-CSR rate and
pressure-gradient functions on the same two-cell synthetic input. This
test does not evolve a Kepler disk.

The previously queued but unstarted v3 fresh `t=0.001→6.5` continuation
was replaced by **437456**, which starts from the passing v4 native pilot
and checks both its one-step audit and v3/v4 physical-field comparison.
The old pending job was canceled at zero elapsed time; existing running
v3 calculations remain untouched. The new job remains an experimental
fixed-mass weighted branch, not a locally conservative ALE/spin solver.
