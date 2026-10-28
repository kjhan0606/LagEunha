# Inner-disk resolution experiment design — 2026-10-27

## Why a new initial condition is needed

The current full-q2 384² IC has 147,456 sites on a 6×6 domain. In the active
disk `0.5 ≤ r < 2`, it has 49,038 sites: 3,931 at `0.5–0.75`, 6,050 at
`0.75–1`, 15,865 at `1–1.5`, and 23,192 at `1.5–2`. The current halo-28
branch changes the periodic domain while retaining these annular counts. The
existing hybrid-polar and ring-spacing generators relocate a fixed total
count; neither tests increased inner resolution. See the
[resolution audit](KEPLER_HALO28_DISK_RESOLUTION_AUDIT_2026-10-27.md).

## Matched local-refinement candidate

Construct a separate W=0 Voronoi IC on the same 6×6 box and halo-28 setup.
Retain the existing exterior coordinates and replace, rather than duplicate,
the sites inside `r=1` with a locally refined layout. A practical first
candidate is half the baseline spacing for `r ≤ 0.75`, with a smooth spacing
transition to the baseline by `r=1`. Leave the geometry at `r ≥ 1` unchanged.
This targets the first active annulus and the low-density center implicated
in earlier close-pair failures, while limiting the added sites.

Area scaling predicts about four times as many sites in the fully refined
region. For the current 147,456-site base, the added count is roughly 20% of
the total after accounting for the transition; the generated IC must report
the exact count and minimum separation. This is an estimate, not a measured
IC result.

The native glass loader accepts either exactly `NX*NY` sites, or fewer sites
when `EUNHA_ALLOW_SPARSE_GLASS=1`; it rejects a glass count greater than
`NX*NY`. Thus a roughly 177,000-site candidate needs a nominal grid with
`NX*NY` above that count (448² provides initial headroom; increase it if the
generated count requires it). Sparse Kepler initialization also requires
`EUNHA_KEPLER_MASS_FROM_CELL_VOLUME=1` and an
`EUNHA_KEPLER_INITIAL_AREAS_FILE` sidecar whose count, box, and area sum match
the refined glass. The accepted [sparse-mass workflow](KEPLER_SPARSE_IC_MASS_CORRECTION_2026-10-27.md)
documents this path. The selected native binary's source provenance must be
checked for these loader features before staging; the full-q2 production
binary is not assumed to contain them. Re-audit the changed `NX/NY` because
they also set the mean-volume reference and cell-list scale.

That sparse loader and periodic-area path already have a native `t=.001`
preflight receipt for 54,312- and 55,931-site glasses. The preserved sparse
glass patch changes only `Exam/KH/util.c` to opt in to `np_file < NX*NY`; the
periodic-area patch also changes only `util.c`, supplying the initial area
before the existing mass rebase. The area-corrected executable used there has
SHA-256
`eff735174de1b487759f168d75e1e317a201692ef8d6cadba03cd03e49afb4b5`.
Those source diffs do not alter the RK4, force, equation-of-state, GPU
tessellation, or timestep code. The active full-q2 reference executable has
SHA-256
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284` and its
run log confirms the variational W=0 pressure branch. The available receipts
support reusing the sparse-area executable for an IC/initialization pilot,
but do not prove bitwise binary equivalence: verify its runtime operator
marker and all physical settings in the pilot log before comparing evolution.

For about 177,000 sites, `NX=NY=448` gives the loader nominal capacity
200,704; the file's count remains the actual particle count. This setting
still changes `meanvol` and the cell-list scale, so verify the post-rebase
initial density and measured runtime rather than treating it as a pure
allocation parameter.

Preserve the same radial density function, pressure, softened gravity, and
orbital velocity. Assign each cell's initial mass from that density function
and its initial Voronoi area, then hold every cell mass fixed during evolution.
Measure the resulting total mass against the base IC; do not silently
renormalize only the refined run, since that changes its physical density
normalization. If a common mass normalization is required, apply it to both
IC definitions and record it explicitly. Keep
all initial weights exactly zero, use no generator steering, and retain the
same native RK4 force operator and 28-cell periodic halo. The baseline and
refined branches must differ only in the inner generator layout and resulting
cell count.

## Cost gate before a long run

The current full-q2 job measures about 2.895 s per RK4 step at
`dt=7.64202e-5`. Halving the smallest cell spacing could nearly halve the
global timestep; adding roughly 20% more sites could raise work per step.
If both effects occur, a fixed-step full-disk run could cost about 2.4 times
the current resolution. This extrapolation is deliberately treated as a
planning bound: the actual timestep and cost must come from a short native
pilot because the smallest cell, local velocity, and mesh geometry determine
the CFL limit.

Before committing to `t=6.5`, run the refined IC only far enough to measure
the accepted timestep distribution, wall time per complete RK4 step, positivity,
and the initial radial-profile checks. Project the full `t=6.5` cost from
those measurements. If the global timestep makes the run several days long,
the next useful development is individual time stepping; do not hide that
cost behind a short density-only result.

## Physical comparison gate

Compare against the full-q2 halo-28 W=0 result at identical times, beginning
with the accepted `t=6.5` checkpoint. Report full density fields and
mass-weighted density error by initial-radius band, especially `0.5–0.75`,
`0.75–1`, and `1–2`; also compare entropy, orbit error, radial motion,
angular momentum, total energy, positive internal energy, exact fixed mass,
zero weights, area tiling, and zero injected energy. The refinement hypothesis
is supported only if inner-disk density/orbit error falls while positivity,
conservation, and the outer-disk errors remain acceptable. The t=6.5 test is
one step toward the requested long-time GIZMO comparison, not its completion.

This experiment isolates resolution. It does not add cell spin or viscosity;
the native spin force/energy closure remains a separate open requirement.
