# Native compression diagnostics — isolated build, not a physical result

2026-10-22, research root `Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics`.
Original source identity remains LagEunha HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`.
The supplied desktop cwd is a different directory; no code there was changed.

## Purpose and current execution

The scientific target remains the historical t6–10 heating/Kepler breakdown and
ultimately the correctly normalized longest GIZMO comparison. The actual new
material trajectory is still accepted step4, t=8.2355880701376864e-5. Neither
this build nor the positive-cubature preprocessing advances physical time.

The original positive-rule producer PID2166274 and its four native ranks were
verified live at 00:31:28 UTC. The processed fraction was 0.3328578255. Dispatcher
PID2351718 was also verified live, waiting for the same complete producer. No
physical transition claim existed. See the exact process start ticks and four
rank progress records in `material_full_compression_diagnostics_v1_ada_mpi4/producer_state_at_build.json`.

The dispatcher and its frozen v2 H200/Ada builds remain the only authorized
downstream execution path for the current paired accuracy/cost bundle. This
diagnostic build is not connected to it and does not release long evolution.

## Implemented in the isolated future build

Directory: `material_full_compression_diagnostics_v1_ada_mpi4`.

The optional native CUDA reductions read the same bounded stage-local cache
after physical pass0 and pass1. They do not evaluate new shapes/fields, change
quadrature, or request an additional force or mass solve. The original physical
19 metrics, MFEvolvedCell/checkpoint ABIs, EOS/stress kernels, mass solver and
paired oracle remain unchanged. The internal, nonserialized cache adds Jgen,
growing from568 to576 bytes/point, or8MiB/rank at its present capacity.
Diagnostic kernel, transfers and output have additional cost that is unmeasured.

For each source cell and each active rule independently, the future output has:

- Its own V0q=sum(w*j0), U0q=P0*V0q/(gamma-1), current Vq, canonical V and U.
- Original canonical V0 inferred from the serialized original initial U;
  this serialization/measure difference is explicitly retained.
- Bulk compression B, signed direct Jensen residual, stable centered Jensen
  remainder, centered lambda variance, centered third absolute moment, and
  centered raw J/Jgen variance. Negative computed gaps/variances are preserved.
- Extrema and original source/compact-positive-point witnesses for lambda,
  J, corrected jhat, Jgen, generating jhat0, and J/Jgen.
- Current/generating delta, mean/volume consistency residuals, sampled Taylor
  lower/upper expressions, and explicitly engineering arithmetic estimates.
- Rule generation, active mode, actual gamma/P0, coefficient/initial-data FNV
  fingerprints, and copied SHA manifests for the generating state, canonical
  input and reduced rule. The ledger time is not claimed to be the RK stage
  physical time; current coefficient fingerprint and stage identify the state.

The small-z remainder uses a binomial series through degree16 with a geometric
tail estimate. The centered gap is recentered for the small residual weighted
mean. This avoids the cancellation in total U/U0 minus a value near1. It does
not turn double Jacobians into exact continuum data, certify libm, measure
irreversible entropy production, or establish the historical heating cause.
Each diagnostic block must visit exactly the physical selected-rule point count.

## Evidence completed

1. Full native SM89 CUDA/MPI/OpenMP build and link: `status.txt` is
   `BUILD_PASS_NO_PHYSICAL_EXECUTION`. CUDA/host compile logs are empty; the
   inherited link has the prior ifx unknown-option warning. Binary/source SHA
   checks and ldd dependency check passed.
2. Static comparison: original full and cached physical EOS/force kernel bodies,
   physical metric ABI, mass code, rule reader and oracle are unchanged.
   `static_operator_comparison.json` and `future_changes.diff` record scope.
3. Host-only arithmetic helper: 293 finite cases compared with210-digit Decimal
   references passed the engineering estimate; one mathematically valid case
   overflows double and four invalid-domain inputs are rejected. The first
   check that treated expected overflow as failure is retained. The final
   `math_validation/result.json` has298 cases and zero classification failures.
   This did not execute the CUDA reduction, CUDA libm, or a physical step.
4. Current producer, frozen v2 Ada/H200 and dispatcher source manifests all
   match: 68+68+68+4 files, zero changes. See `live_frozen_source_integrity.json`.

## Conditional use and remaining work

Do not launch this build before the current paired accuracy/cost result and
the ensuing decision. No diagnostic physical execution was performed here.
The runtime feature requires the existing continuous cached path, defaults off,
and is capped at16 evaluations per rank. It needs:

    GFS_MATERIAL_COMPRESSION_DIAGNOSTICS=1
    GFS_MATERIAL_COMPRESSION_CANONICAL_INPUT=<original initial-ledger run>
    GFS_MATERIAL_COMPRESSION_CANONICAL_SHA256=<this directory>/canonical_input.sha256
    GFS_MATERIAL_COMPRESSION_RULE_SHA256=<SHA manifest for exact loaded rule files>
    GFS_MATERIAL_COMPRESSION_MAX_EVALUATIONS=<1..16>

An eventual launcher must externally verify the source, binary, immutable,
canonical and exact rule SHA manifests before native execution. The native
diagnostic copies/binds those manifests and checks identities/counts; it does
not itself implement SHA256. The old original U0 and the quadrature transition
ledger are never overwritten by its rule-specific diagnostic U0q.

GPU numerical agreement, actual overhead, continuous enriched-map/local-degree
guard, spatial representation lifetime, W0 reduction, and t6/8/10 plus longest
GIZMO validation remain incomplete. The last accepted physical state is retained.
