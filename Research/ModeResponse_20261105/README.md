# Finite mode and known-wave response audit

See ../../Docs/UNIFIED_HYDRO_MODE_RESPONSE_2026-11-05.md for mathematical
contracts, results, and limits. Original native 1D/2D/3D geometry and the saved
N128 coefficient reference are reused. Python only evaluates finite matrices,
known Fourier probes, and bounded spectral diagnostics. No fluid is evolved.

Authoritative results are material_volume_modes_v2/results. Independent
verification is material_volume_modes_v2/verification/result.json. The failed
v1 analysis is retained. The v1 mass-array adapter shape error preceded any
native spectral result. An explicit unequal-mass two-cell check validates the
corrected two-component mass factor before v2 runs. No coefficient or PDE
simulation was repeated.

C passes sampled long-wave and transverse consistency gates but worsens
pressure force on the saved deformed geometry. The native positive acoustic
frequency is 49.56 times the B frequency. Independent reassembly and bounded
power iteration reproduce the candidate value. This is neither a growth rate
nor a full KDK/background stability qualification. The dominant candidate mode
is localized around three cells with volumes near the mean. All spectral
verification work remains inside the combined preregistered caps.

Compact source snapshots, outputs, commands, validations and SHA256 evidence
are archived. Large inputs remain on GPFS with hashes. Existing CPU coefficient
cost/tolerance failures persist. No solver tuning, limiter or candidate PDE is
released. The full physical accuracy, dimensional PDE suite and matched GIZMO
throughput objectives remain open.
