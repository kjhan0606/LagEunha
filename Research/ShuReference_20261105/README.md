# Original Shu--Osher IC repair and independent Euler reference

See `../../Docs/UNIFIED_HYDRO_SHU_REFERENCE_2026-11-05.md` for the result,
limits and exact comparison scope. The original common initializer used
sin(5x), while the original legacy problem uses sin(5*pi*x). No previous
accepted state is rewritten. The native common initializer is corrected;
the operator and integration method remain unchanged.

`native_run/result.json` is the final same-IC comparison to Athena N12800.
`native_run/result_reference_N6400.json` preserves the intermediate reference
level. `initial_condition_mismatch_raw.json` is explicitly excluded from
accuracy claims by the adjacent explanation JSON. Reference uncertainty is
reported, not set to zero. The density plot uses actual FV data.

The 25 native PDE processes are terminal: four Athena++ references through
t=1.8, five corrected common implicit-KDK Shu resolutions through t=1.8,
and parent/current controls for the other eight N200 problems at their
original endpoints. The latter eight final data pairs are byte-identical.
`publication_validation.json` records endpoint/source checks. Historical
`running.json` files are launch provenance, not current live handles.

The external reference is pristine upstream Athena++ commit
`37869da22cd357c2fbdd056fbb3781cdd4f430da`. No uncommitted changes from the
user's existing Athena checkout enter it. Its BSD license is preserved as
`LICENSE.athena`; the upstream problem/EOS patches, originals, configured
Makefile/defs and source hashes are included. The reference has native
C++ HLLC/PPM/RK3 dynamics. Python only schedules and integrates saved norms.

Most raw tables, source snapshots and native logs are losslessly gzipped,
preserving all original bytes, including native whitespace. Every archive
item has compressed and original SHA256 metadata in `manifest.json`.
Executables remain under the GPFS run root with binary SHA256s recorded.
The entire external source tree and object files are not duplicated here.

The authoritative run directory is
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/shu_independent_athena_v1`.
With that preserved directory, analysis/figures can be replayed without a PDE:

```bash
python3 /gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Exam/HydroCore/shu_independent_reference.py \
  /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/shu_independent_athena_v1 --plot
```

For a fresh rebuild, use a new GPFS directory, pristine upstream commit,
the two recorded Athena patches and configure/build commands. The common
native frozen dependencies and exact build argv are included separately.
`--run` creates the registered reference cases and `--run-native` creates
the native suite, refusing existing case directories. The current replay
source is `shu_independent_reference_v3.py`; earlier driver versions are
preserved to explain the input mismatch and the one extra reference level.

This is not a result for the newly planned FV--ALE candidate. Common
multidimensional physics and matched GIZMO cost remain open.
