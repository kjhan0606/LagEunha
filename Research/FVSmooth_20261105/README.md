# FV smooth-flow evidence and retained failure

See `Docs/UNIFIED_HYDRO_FV_SMOOTH_2026-11-05.md` for the full interpretation.
These are actual native C/OpenMP trajectories and read-only native
actual-polygon references. Python schedules and assesses them.

The registered HLLC/HLL candidate completes Gresho64/128 and its 64
half-step control, but fails the cold-pressure criterion. The earlier
failed endpoints and later same-time controls remain distinct. Native
thermal chain-rule audits and restart replay evidence are included.

`kep_flux_v1` contains only the initial optimized C/C++ algebra of the
replacement core. Its later actual coupling is archived separately in
`Research/KEPNative_20261105`. Do not treat this algebra archive as a
complete PDE or GPU qualification.

`manifest.json` gives source/stored SHA256 values and compression flags.
Gzip files preserve the raw bytes, including source whitespace. Native
executables are excluded; source/build manifests identify their hashes.
Historical `running.json` launch records are not live job status.
Original binaries and uncompressed data remain under
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`.

The `gizmo_source_receipt.json` records inspection of retained pristine
public GIZMO source. No new matched GIZMO run or speed ratio exists here.
