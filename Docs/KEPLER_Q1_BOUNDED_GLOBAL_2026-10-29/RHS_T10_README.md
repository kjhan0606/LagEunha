# Native t10 pressure-force and thermal-work snapshot

The exact accepted t10 W0/fixed-mass checkpoint was read in a separate diagnostic-only native v2 CUDA+MPI run on Ada. The program computed K1 for all 90,352 owned sites, saved one 52 MB diagnostic state, emitted eight selected-ID TSV rows and exited before any time evolution. The selector was fail-closed: one exact step, one RHS stage, eight unique IDs, no face/tessellation dumps and a file-size ceiling. `RHS_T10_PROBE_INPUT.json` specifies the choices, deliberately including core/adjacent extremes and two disk sites; they are not a representative random sample.

For each of those eight sites, independent complete bounded Voronoi polygons reproduce the native pressure acceleration and pressure-work RHS within the preregistered absolute tolerances; `RHS_T10_SELECTED_AUDIT.json` includes the error and bounds on omitted roundoff-length edges. The selected cells alone do not establish population behavior.

For the full original-ID core plus disk (`r0<2`, 51,926 cells), a separate independent unweighted Voronoi reconstruction computes all touching finite faces and reconstructs native variational pressure force and adiabatic thermal work. Native cell areas agree within 5.76e-8 maximum relative. The independently reconstructed force and thermal work match the eight native sample rows to maxima 9.39e-15 and 6.52e-15 respectively. `RHS_T10_POPULATION.json` records all radial bands, tails and source hashes. The population reconstruction ran on Ada in 2.24 seconds with Python 3.13.11, NumPy 2.3.3 and SciPy 1.18.0.

| Original region | Cells | Cell median pressure acceleration / central gravity | Mass-mean abs(Udot)/(U Omega) | Cell p99 abs(Udot)/(U Omega) |
| --- | ---: | ---: | ---: | ---: |
| r0<.1 | 134 | 3.0402e-5 | 0.23598 | 0.88917 |
| .1–.5 | 3,112 | 7.9513e-8 | 0.012659 | 0.39965 |
| Disk .5–2 | 48,680 | 5.3760e-8 | 0.0034907 | 0.015311 |

Omega is the local Plummer circular frequency at the current site position. The pressure/gravity ratio is ill-conditioned near the central generator where gravity approaches zero; its core maximum is not representative. Absolute pressure accelerations are also reported. The thermal ratio is an instantaneous compressional/expansional work rate, **not entropy production or a finite-interval heating measurement**. This is one already core-damaged state, not a core-accurate solution or proof that either component remains slow throughout a macrostep. It does show that small pressure acceleration alone is insufficient to assume thermal work is slow: a proposed gravity subcycling scheme must resolve and converge pressure/geometry/work coupling independently.

A first analysis used `np.linalg.norm` rather than the native gate's canonical `np.hypot` for original radii. Five IDs on exact radial bin boundaries were consequently assigned differently. That result and its exact source are preserved as `RHS_T10_NONCANONICAL_INVALID.json` and `rhs_t10_noncanonical_source.py.txt`; `RHS_T10_NONCANONICAL_RECEIPT.json` explains the correction. The table and current script use the canonical gate and reproduce all 48,680 disk IDs. Only read-only analysis was repeated; no native evolution or probe was rerun.

`RHS_T10_NATIVE.log` is a display copy with trailing spaces removed. The unchanged native log and its SHA256 are recorded in `RHS_T10_SELECTED_AUDIT.json` under the original GPFS path.
