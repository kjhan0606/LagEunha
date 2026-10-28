# TCF pair exchange and thermal identity audit

Source inspected: Exam/exam.c, active TCF assembly near lines 9122-9140.
Using 221159 t6 recorded state/face rows, with explicitly synthetic tangential
offsets (actual face centroids are absent from the dump), the full pressure plus
recorded-Q expressions satisfy owner reversal exactly for force, energy,
volume and heat. A common velocity boost (7,-3) satisfies energy covariance
to 2.17e-19 and volume covariance to 1.67e-16 absolute error.

The per-face thermal identity is

    Edot_i - v_i.F_i + P_i Vdot_i - H_i = A P_i v_i.e

and holds to 1.73e-20 in this probe. The right-hand side cancels only after
summing a closed cell. Thus checking zero per-face residual would incorrectly
reject the implementation. Geometry closure must be checked separately.

CORRECTION: the initial radius filter mistakenly used (0,0), whereas the dump
center is (3,3), as used in existing pair-force and face-heat analyses. The old
JSON is retained as superseded evidence; its radial classification is invalid.
With the correct center, r<1.9 contains 11498 cells, max |P v.sumAe|=9.36e-21
and RMS 6.20e-22. Max norms of sumAn and sumAe are 1.28e-16 and 1.24e-16.
The global maximum 8.14e-7 belongs to cell 36858 at r=4.09375. All ten largest
terms are at r=4.0625 or 4.09375. This audit provides no support for interior
face-closure error as the instability mechanism. These are diagnostic rates,
not accumulated energy errors or a complete boundary-condition audit.

This is an independent algebra implementation, not production-source replay
or CUDA/MPI validation. Actual upstream geometry, QBP, HG, RK4 and physical
accuracy remain outside this gate. Existing production runs are unchanged.

Corrected results: canonical_tcf_pair_identities_corrected_2026-10-18.json.
Superseded results: canonical_tcf_pair_identities_2026-10-18.json.
Script: audit_tcf_pair_identities.py. Remote root:
/gpfs/kjhan/LagForce/Codex_tcf_pair_identities_corrected_20261018.

## Follow-up across saved cases

`audit_tcf_closure_history.py` and
`canonical_tcf_closure_history_2026-10-18.json` check t2/t4/t6/t10-named dumps
with the correct center (case labels, snapshot times not independently read).
Across these cases, interior sumAn/sumAe maximum norms remain below 1.31e-16.
Inner-region thermal closure RMS ranges from 6.40e-22 to 6.70e-21 while
recorded TCF heat RMS grows from 2.04e-10 to 1.07e-7. The closure contribution
is too small to explain that measured heating. Exterior maxima are unchanged
across all four saved cases, consistent with fixed exterior data but not by
itself proof of the boundary implementation.

This eliminates the newly suspected interior closure defect in these dumps.
The established distinction between geometric cell-volume evolution and
material-volume evolution remains. Next backend work should test original
TCF assembly with actual face centroids and upstream QBP inputs, rather than
change correct geometric normals or remove physical pressure work.
