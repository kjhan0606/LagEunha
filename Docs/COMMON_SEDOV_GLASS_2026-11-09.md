# Sedov glass-layout diagnostic — 2026-11-09

The user proposed glass ICs instead of a regular grid. This is an initial spatial-sampling control for the unresolved2D Sedov morphology. Preserve all completed square/disk runs and common hydro; do not start3D or add a force/viscosity correction. The current disk128 evolution has finished and is undergoing CPU polygon analysis; this new layout diagnostic uses no GPU or evolution.

## Existing implementation audit

`Exam/KH/generate_glass_2d.py` provides periodic Lloyd relaxation with3×3 Voronoi tiling. Its default initial points are a regular grid plus independent displacements in[−.25h,+.25h], despite a nearby comment describing50% jitter. `Exam/Sedov2D/util.c` has a separate glass option starting from a grid with default±.15h jitter and a pinned central anchor; its comment claiming a uniform-random start does not match that branch. The legacy Sedov implementation DOES recompute m=rho*V and U=P*V/(gamma−1) after relaxation (around line810); do not claim that it simply retains mean-volume masses. The native common-method comparison should reuse the consistent actual-cell initialization principle, not silently inherit the old anchor or source deposition.

## Completed layout-only control

`common_sedov_glass_20261109_v1/layout_audit.py` imports the existing KH periodic Lloyd function unchanged, compares the original KH jitter seed to a uniform-random seed, and applies100 iterations at4096points. Seed42, same periodic unit box and relaxation routine. Original source hash is in plan.json; positions and actual periodic Voronoi areas/centroids are saved in three NPZ files. No new glass algorithm or hydro method is introduced.

For unique periodic neighbor pairs separated by less than1.3h (h=1/sqrt(N)), define A_m=|mean exp(i*m*theta)|. This measures initial bond orientation, not a density Fourier spectrum, shock contour amplitude or physical error.

| Layout | A4 | A6 | Cell-area CV |
|---|---:|---:|---:|
| Regular square | 1 | ~0 | 0 |
| Existing KH grid-jitter + Lloyd | .0670190 | .0295404 | .0482159 |
| Uniform random + same Lloyd | .00434157 | .0384983 | .0698449 |

All periodic cells have positive area, sum1 within1e−11. The random-start layout has much less global fourfold bond alignment in this seed, but finite local/order fluctuations remain; Lloyd is not asserted equivalent to a gravitational glass and one seed is not an ensemble result. Greater cell-area scatter must not be converted into a density perturbation by assigning identical masses while assuming equal volumes.

## Next implementation contract

Use the random-start relaxed layout as a controlled alternative to the existing regular disk IC. With W=0 and rho0=1, assign each fixed material mass from its actual cell area, m_i=rho0 V_i. Set background thermal energy from P0 V_i/(gamma−1), and integrate the same finite disk source into each actual convex polygon. Then evolve those masses unchanged using the current native C/CUDA solver. Before evolution, independently certify cell geometry, total mass/excess energy, uniform background pressure/density and zero velocity; validate source integration against an independent calculation. Keep the same source radius, target time, EOS and solver; source shape and layout must not be changed together. Multiple predeclared seeds are needed before a general isotropy claim.

At this report only point layouts and geometry diagnostics are complete. No glass native checkpoint has been written and no glass hydro evolution has started. No old result is replaced or relabelled. Artifacts: plan.json, layout_status.json, layout_report.json, regular.npz, legacy_jitter_lloyd.npz, random_lloyd.npz. No commit/push.


## Actual-cell native IC completed and queued

The seed42 native64² checkpoint is now prepared. Existing periodic power-cell reconstruction and independent SciPy periodic Voronoi areas agree within3.12e−14 relative. Source deposition uses exact circle/edge sector-triangle integrals and is checked independently by adaptive vertical-chord quadrature; largest intersection-area difference1.0842e−19. Four mathematically empty overlaps had negative geometric roundoff down to−1.0842e−19, within64*machine-epsilon*cell-area; these are explicitly recorded and set to their nonnegative zero value. No energy renormalization is applied.

Total volume and mass are1, excess explosion energy.9999999999999989. Each initial mass equals its actual cell volume, giving rho=1 exactly. Background pressure error≤2.12e−22 absolute. Native codec round-trip retains coordinates/masses/volumes exactly, U agrees with the prepared value exactly, and all velocity/weight momenta remain zero. Initial masses vary from cell to cell but remain individually fixed thereafter. Table row count, truncation, extra trailing data and nonconservative energy are rejected before writing a snapshot.

An initial attempt to add a table branch inside the regular IC loop changed KH vy by up to1.73e−18 due to changed compiled arithmetic and failed strict byte-parity. No evolution used that candidate. The tested regular writer was restored byte-for-byte (SHA25210e3f807f19cc716f53ecda3825f91337c44c3491053a762c6fdc7da92f50). The final implementation is the small `Exam/NativeU64/sedov_cell_initial.c` adapter using the existing native snapshot codec and requiring a pristine t=0 unit-box Sedov template. Six default IC checkpoints and the regular disk checkpoint reproduce the originals bitwise. Failed candidate and comparisons remain archived; no earlier hydro result is changed.

Artifacts: `initial64/{initial.dat,state_ids.bin,cells.txt,geometry.npz,audit.json}`, `ic_adapter_v1/prepare_ic.py`, compiled adapter and build commands, `table_adapter_decision.json`. Initial checkpoint SHA b5e3a7df3c116ff33a8b761c976987281161f2725708f595e13ba9c101d7661b; adapter source SHA ce36a7b160ce90849140a6863cd539dd39e07ffc5edd3b973232a26e0f8e8995.

Detached supervisor487245 has queued one64² glass evolution through t=.1 after current Gresho128 completes and its actual final native process exits. It uses the same11a0234… hydro binary, C=.4, dtcap.000125, endpoint tolerance1e−11, disk radius and periodic boundaries as the regular64 disk. A finite six-hour dependency wait checks process identity, not just status files. No live calculation is stopped or restarted. Checkpoints every.05, maximum4 including input (~11.4MB);149.09TB free at preparation. Effective params will be `run64/params.dat`. Final geometry/energy gates and the same contour/profile analysis follow; failure stops the comparison without a patch or reset.

At queue registration, glass physics has not yet started. This is a single-seed comparison; it cannot establish ensemble isotropy or complete2D qualification.


## Finite-cell source projection recorded before evolution

Degree4 polynomial integration over the actual cell polygons gives the piecewise-constant deposited source centroid(.4999529314,.5001049911), offset.007364 mean cell spacings from the continuous center. Its second radial moment is10.7600% above the continuous disk value (regular64 disk was12.264% above). Normalized quadrupole components H2/R2=(−.0109058,+.00063736);28cells contain more than1e−12 source energy. These are finite-cell projection effects at this one seed, despite the exact geometric disk-overlap energies. They are recorded in `initial64/projection_moments.json`; no recentering, renormalization, seed selection or compensating force is applied. A reduced global bond A4 does not guarantee an isotropic projected explosion or a lower hydro error.


## Completed glass evolution and inspected morphology

The seed42 random-start Lloyd calculation reached t=0.1 with the unchanged common C/CUDA binary and original planned settings. 844 accepted steps,16 retries; all material masses are bitwise fixed, U/V remain positive, total area1, energy drift−3.33923e−13 relative and maximum accepted volume mismatch9.90770e−12. Solver wall time215.776s was measured on shared Ada and is not a matched GIZMO speed certificate.

| Layout (4096 cells) | Density L1 vs point source | Pressure L1 vs point source | rho=2 contour radius RMS/mean | Fourth angular amplitude |
|---|---:|---:|---:|---:|
| regular_disk | 10.936165% | 9.577644% | 0.903961% | 0.135684% |
| glass | 7.209467% | 4.450385% | 0.763528% | 0.202310% |

Actual-cell density panels were visually inspected: the repeated diagonal rays and regularly arranged low-density pockets of the regular grid weaken in the random-start result. Irregular low-density patches and a cell-scale rough shock ring remain. The rho=2 radius RMS improves, but the fourth-mode amplitude does not; at rho=2.5 the radius RMS is essentially unchanged (.768707% regular, .769390% glass). Do not generalize one threshold, one seed or the initial bond A4 as proof of isotropy. Source projection and fixed material masses differ between layouts even though the continuous disk source and initial density/pressure are matched. No physics terms, viscosity, limiter or tunable coefficient changed.

The analytic reference is an ideal point explosion, whereas BOTH compared ICs use a finite DISK. The inherited assessment JSON caveat saying 'finite square source' is a stale label; it does not describe these runs. The radius/profile calculation is shape-independent point-source postprocessing. Preserve the raw reports, and use this explicit source provenance when interpreting them.

The first density render omitted periodic translated polygon pieces outside the displayed unit box, leaving white strips at its edges. This was a plotting omission, not a vacuum or missing hydro cell. `gresho_glass_evidence_20261109_v1/periodic_figures.py` now clips translated saved polygons to the displayed box and checks their summed area equals1. The corrected `glass_figures_periodic_v2/density.png` and PDF were visually inspected. They use unsmoothed cell densities, identical0..5 colors, no out-of-range cells, and the same dashed point-source circle. Original plots remain preserved. The corrected image is authoritative for display.

Glass improves the radial/profile comparison in this controlled seed but does not finish Sedov qualification or remove all morphology error. A predeclared seed/resolution comparison is needed before adopting it as a general remedy; no3D run is justified from this result alone.


## Predeclared refinement and seed comparison launched

`common_sedov_glass_refinement_20261109_v1/plan.json` fixes the next cases BEFORE their results:128² seed42,64² seed7,64² seed19. All use100 iterations of the same existing periodic Lloyd routine, actual-cell initial masses and disk-source integrals, the same common native binary, C=.4, source width.0625 and target.1. The timestep cap scales as.000125*64/N, matching existing regular-disk controls. Seed42 at64 and regular64/128 are reused from completed results, never rerun or discarded. This is a small exploratory seed set, not statistical proof of isotropy.

The batch first independently checks actual power areas against SciPy periodic Voronoi geometry and source-circle intersections against vertical quadrature. It uses the previously verified frozen native cell-table codec; inherited codec parity is explicitly distinguished from newly checked IC geometry and round-trip fields. No source recentering, energy renormalization, best-seed selection, viscosity change, or problem-specific force is applied.

Supervisor558352 was launched after the8-case acoustic batch completed. It prepares128 ICs on CPU22,23, then runs each hydro case sequentially in the released third Ada slot on CPU12,13. Existing RT/cylinder are untouched. The12-checkpoint bound includes copied inputs; estimated total68.42MB,149.18TB free at launch. Maximum temporary dense geometry4.30GB, available host RAM816.85GB. Effective per-case params/env and child handles are recorded by the reused runner. Any invariant failure stops the batch for inspection, without automatic retuning or restart. Raw cases and scripts are preserved.


##128 glass endpoint and inspected density image

The predeclared seed42 run at128² completed t=.1:1708 accepted steps,30 retries,888.339s on shared Ada, fixed masses, positive U/V, total volume1, energy drift3.591e−13 and maximum volume mismatch9.709e−12. The subsequent64² seed7 run started automatically; seed19 remains next in the same finite sequential batch.

At128² the point-source density L1 is4.21033% for glass versus6.39036% regular; pressure L1 is2.51401% versus4.28370%. The rho=2 contour radius RMS is.393316% versus.398678%, essentially similar, while the fourth angular amplitude is.0573231% versus.129706%. Thus profile and directional-pattern improvements do not mean a large improvement in every shape statistic.

`common_sedov_glass_refinement_20261109_v1/figures128/density.png` and PDF were rendered from actual saved power polygons, with periodic translated pieces, unit area, identical0..5 density colors and the same point-source reference circle. Visual inspection confirms the regular grid's diagonal rays are much weaker in glass, while irregular interior patches and shell noise remain. This is still one128 seed, and the reference is not an exact finite-disk solution. The other predeclared seeds must be retained regardless of whether they improve or worsen the result.


## Completed predeclared 64-square seed comparison

All predeclared additional runs completed at t=0.1:128² seed42 and64² seeds7,19, supplementing the original64² seed42. Supervisors and native processes for this batch are terminal; no restart or seed selection was performed. All retain fixed cell masses, positive states and total energy residuals below7e-13. The64² comparison table and four-panel actual-cell density figure are in `common_sedov_glass_refinement_20261109_v1/ensemble64/report.json`, `density.png`, `density.pdf`; the plotting script and input hashes are saved. Figure visually inspected: no missing periodic edge pieces, shared color limits0..5, no smoothing; the three glass layouts reduce coherent grid-aligned rays but retain irregular density structure.

|64² layout|Density L1|Pressure L1|
|---|---:|---:|
|Regular grid|0.10936165|0.09577644|
|Glass seed42|0.07209467|0.04450385|
|Glass seed7|0.06837968|0.03942464|
|Glass seed19|0.07646814|0.04768094|

All three glass seeds improve these profile norms against the SAME point-explosion reference, which is not the exact finite-radius disk IC solution. Density error is30.1–37.5% below regular, pressure50.2–58.8% below regular. This three-seed sensitivity check is not a statistical confidence interval, and it does not establish improvement of every angular mode or density threshold. At128² seed42 the density/pressure norms are0.04210333/0.02514015 vs regular0.06390364/0.04283703, but rho=2 contour angular RMS remains nearly unchanged (0.00393316 vs0.00398678).

Here “glass” means random-start Lloyd-relaxed periodic generators,100 iterations, not gravitational glass. Uniform initial density is represented by m_i=V_i, so initial cell masses are unequal and then remain fixed. This is an IC-layout change with unchanged hydrodynamics; it is evidence for reducing lattice imprint, not a replacement for operator consistency checks or proof of convergence in3D.
