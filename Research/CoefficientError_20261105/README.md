# Frozen coefficient error and spatial response

See ../../Docs/UNIFIED_HYDRO_COEFFICIENT_ERROR_2026-11-05.md for identities,
results, interpretation, and limits. The original native C/OpenMP probe made
one coefficient capture under the unchanged cap, tolerance, and shared budget.
Both fields failed convergence. No native PDE steps were performed.

The authoritative analysis is analysis_v2/propagation.json. The first analysis
is retained. The revision uses the same symmetric analysis edges throughout
the commutator and bounds and adds the exact/deviation response decomposition.
The native capture and existing global reference were not rerun. NumPy/SciPy
only analyze frozen operators and independent small-matrix mode responses.

Tree residual-flow estimates do not use Ritz estimates as certified spectral
gaps and are not interval-arithmetic proofs. Reference uncertainty is explicit.
No new tolerance, successful cost qualification, or candidate PDE is released.
Accurate coefficients improve the exact Gresho velocity response but amplify
the deviations in the saved original-operator state. The latter response is
not proof of instability of a candidate trajectory.

Source snapshots, commands, output, and SHA256 provenance are archived.
Large binary/array artifacts stay on GPFS with hashes in manifest.json.
Production physics and the manuscript benchmark remain unchanged. The full
dimensional accuracy, remaining 2D suite, and matched GIZMO speed goals remain open.
