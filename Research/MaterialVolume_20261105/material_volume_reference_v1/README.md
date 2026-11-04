# Independent frozen-graph coefficient reference

`Exam/HydroCore/material_volume_reference_audit.py` reads the direct native
128² Gresho t=0.5 coefficient export. A pinned sparse factorization supplies
arithmetic-mean-zero reference coefficients. Ritz values include measured
residuals and are explicitly not certified spectral bounds.

The archived result is from a completed reference calculation. The saved
matrix, coefficient, and RHS arrays were checked again before publication.
`manifest.json` records the original input and output paths, sizes, and SHA256
hashes. Large arrays remain on GPFS. The manifest captures the source after the
run and does not claim a pre-run frozen source snapshot.

No fluid evolution, production solver change, new tolerance, or successful
candidate cost qualification follows from the reference result. The next
residual-to-physical-error analysis remains pending. The current bounded CPU
coefficient implementation remains held.

To reproduce, invoke the script with the recorded input and a new output
directory. The script refuses to reuse an existing output directory. Timing
and Ritz values can vary between runs. Reference timing is not production
performance evidence.
