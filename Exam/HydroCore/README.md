# Common hydrodynamic operator audits

The shared headers are experimental. Algebraic identities, successful 1D
evolution, and native saved-state RHS checks are separate validation scopes.
See `../../Docs/UNIFIED_HYDRO_PROGRESS_2026-11-04.md` and
`../../Research/UnifiedHydro_20261104/README.md` for the current limitations.

Run from a fresh directory under `/gpfs`. Set `repo` to the absolute checkout
path and `run` to a new calculation directory. Never use `/scratch`.

```bash
repo=/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1
run=$(mktemp -d /gpfs/kjhan/LagForce/common-hydro-audit.XXXXXX)
cc -O3 -std=gnu11 -Wall -Wextra "$repo/Exam/HydroCore/pair_audit.c" -lm -o "$run/pair_audit"
cc -O3 -std=gnu11 -Wall -Wextra "$repo/Exam/HydroCore/fit_audit.c" -lm -o "$run/fit_audit"
"$run/pair_audit"
"$run/fit_audit"
cc -O3 -std=gnu11 -Wall -Wextra "$repo/Hydro1DExam/laguerre_sod.c" -lm -o "$run/hydro1d"
mkdir "$run/blast" "$run/shuosher"
cp "$repo/Hydro1DExam/blast/ref.dat" "$run/blast/ref.dat"
cp "$repo/Hydro1DExam/shuosher/ref.dat" "$run/shuosher/ref.dat"
cd "$run"
./hydro1d unified_implicit 200
```

The 1D command runs stationary contact, moving contact, Sod, blast, Noh,
Lax, double rarefaction, collision, and Shu--Osher with one coefficient set.
The `unified` argument selects the explicit dissipation candidate.
Use separate run directories to preserve outputs from each selection.
The historical Shu--Osher reference has an upstream-boundary defect. Norms
against the historical reference do not establish accuracy. Use the spatial
self-convergence audit recorded in `shu_reference_audit.json` and obtain an
independent reference before claiming verification.

`geometry_audit.c` uses the project's actual Voro implementation.
`device_audit.cu` evaluates shared functions on a CUDA device.
`integrator_probe.c` replays one state with the RHS compiled from the current
headers. Reproducing a rejected historical candidate requires the matching
historical source snapshot in the GPFS experiment directory, rather than the
current headers.

`implicit_probe.c` is a frozen-state 2D sparse linear-algebra experiment.
The probe is not a hydrodynamic timestep and does not solve endpoint active
signs. The binary input layout is explicitly declared by the source structs.
The archived native build invocations and results are in
`../../Research/UnifiedHydro_20261104/results`.

The native CUDA module is `../NativeU64/gfs_unified_gpu.cu`. The adapter
supports one-rank W=0 fixed-wall or explicitly selected periodic domains.
Periodic images are aliases of owned cells rather than independent unknowns.
Enabling `GFS_UNIFIED_HYDRO=1` alone does not qualify an arbitrary boundary,
MPI decomposition, or 3D evolution.


Later bounded native coupling is recorded in
`../../Docs/UNIFIED_HYDRO_NATIVE_IMPLICIT_GATE_2026-11-04.md`.
The persistent Kepler heating remains unresolved. Subsequent actual periodic
2D contact, steady shear, and planar shock runs are documented in
`../../Docs/UNIFIED_HYDRO_PERIODIC_2D_2026-11-04.md` with compact evidence in
`../../Research/Periodic2D_20261104`. Axis-aligned Sod has unresolved transverse
velocity growth. Full 2D benchmarks, actual 3D evolution, and matched GIZMO
performance remain incomplete.

The actual 1D stage and time-refinement checks are `implicit1d_audit.c`.
`implicit_graph_audit.c` compares the common graph solve against the 1D
tridiagonal solve and checks conservation in dimensions 2 and 3. Link the
latter with `../gfs_viscosity_implicit.c` and OpenMP.
`periodic_viscosity_audit.c` checks the same solve across a periodic cut,
including angular-momentum flux, in dimensions 2 and 3.
`viscosity_failure_replay.c` reads a captured half-step graph input. The input
is a diagnostic record, not an accepted hydrodynamic restart.
`exact_sod_sample.c` samples the analytic Sod solution for point-based norms.
`closure_channel_audit.c` reads the frozen native binary format from
`implicit_probe.c`. Named circular-gradient and area-fit alternatives are
causal diagnostics only. `scale_separation_audit.c` links the actual
`../../Voro/voro.c` with XYZDBL and GOTPM definitions and checks fixed smooth
and checkerboard velocity fields at three spatial resolutions.

The source and compact data inventory for the later stage is in
`../../Research/NativeImplicitGate_20261104`.
