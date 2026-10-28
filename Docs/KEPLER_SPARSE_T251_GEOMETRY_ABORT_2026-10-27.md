# Sparse q2 Kepler: paired GPU clipping abort near t=2.51

## Observed endpoint and retained state

The accepted CUDA+MPI+OpenMP A10 continuations **437253** (wide, 55,931 sites)
and **437254** (narrow, 54,312 sites) both ended `FAILED 1:0` before their
planned `t=10` endpoint. Wide failed while building RK4 stage geometry from
accepted step-start `t=2.5176961142503815`, and narrow failed in the final
volume tessellation just after `t=2.50986`. Both logs report
`[GPU_CLIP_CAPACITY] nc=48 capacity=48`, followed by a geometry abort.
Neither is a measured disk physics failure or a `t>6` result. The unchanged
accepted wide executable SHA-256 was
`eff735174de1b487759f168d75e1e317a201692ef8d6cadba03cd03e49afb4b5`.

Wide's complete accepted failure-start snapshot was saved as
`khdiag.step001329.start_geometry.dat`. It is explicitly the **accepted
step-start**, not the failed RK stage. An independent native-versus-exact-orbit
[audit](KEPLER_SPARSE_WIDE_NATIVE_BALLISTIC_T2517696_ACCEPTED_ABORT_2026-10-27.json)
checks complete IDs, positive state, fixed mass, W=0, periodic area 36, and
`rho=m/V`. It finds original-disk mean absolute native density change
`1.39326e-4`, exact-orbit density change `1.39420e-4`, and pointwise
native-versus-orbit residual `6.68859e-7` at this retained state.

The latest regular outputs are wide `t=2.4761277758878566` and narrow
`t=2.4753977144699952`. Their independent [wide]
(KEPLER_SPARSE_WIDE_NATIVE_BALLISTIC_T2476_A10_2026-10-27.json) and
[narrow](KEPLER_SPARSE_NARROW_NATIVE_BALLISTIC_T2476_A10_2026-10-27.json)
receipts give disk native/ballistic mean absolute density changes of
`1.43444e-4`/`1.43519e-4` and `1.43986e-4`/`1.43525e-4`, respectively.
Pointwise residuals are `6.50818e-7` and `1.47558e-6`. Thus both disk
fields remained close to exact-orbit geometry **before** the GPU abort.

The [wide](KEPLER_SPARSE_WIDE_SEAM_PROVENANCE_T2476_A10_2026-10-27.json)
and [narrow](KEPLER_SPARSE_NARROW_SEAM_PROVENANCE_T2476_A10_2026-10-27.json)
cell-ID provenance receipts show a new limitation relative to the earlier
`t=1.65` check: **24 initially outer/seam-near cells, fixed mass
`0.00584865`, have moved inside the `r=2.35` guard** by `t≈2.476`.
None has a generator center inside the original `r<2` disk yet, and the
minimum center gap to its outermost original-disk generator is about `0.256`.
Seam-tagged mean absolute density change is `1.107` (wide) and `0.987`
(narrow). Center crossings alone do not measure face-mediated pressure,
work, or torque transfer; a late disk isolation claim remains unproven.

## Corner-cap diagnosis and isolated probe

The native GPU raw-face clipping loop starts with a box and processes
halfplanes in linked-cell order. It merges consecutive collinear segments
with the same support **only after** all clipping. A read-only, independent
SciPy periodic Voronoi reconstruction of the accepted wide failure-start
positions and narrow latest positions gives maximum final cell vertex
counts **10** and **11**, respectively, with every primary site bounded.
The reconstruction included 55,692 wide and 53,004 narrow periodic images
selected within 2 units of the relevant box faces. The true final cell
complexity is therefore far below 48; the GPU cap is encountered in an
intermediate clipping polygon. That observation suggests support-fragment
accumulation, but a successful native replay is still needed to verify the
fix and detect any other clipping error.

An isolated copy of the accepted source at
`/tmp/Kepler_sparse_clipcompact_20261027` now applies the **same existing
support-identity coalescence rule during exact clipping when an intermediate
polygon reaches 24 vertices**. No physics, mass, weight, timestep, halo,
or corner capacity was changed. The exact [source patch](KEPLER_SPARSE_CLIPCOMPACT_2026-10-27.patch)
is relative to the preserved accepted source overlay. The only edit is in
`Exam/exam_gpu.cu`, SHA-256
`f00ee518d0ca02ecb68eb88c22efd1b974d2eb30611475c12d97516ef14a33cd`.
The new native executable is
`886638100a4ae057ccd6c3624c1b712141dea0656ccaa791934d0f00d44834fc`,
with `sm_86` CUDA cubin, GPU tessellation symbol, and no missing libraries.

Short A10 job **437277** restarted exactly once from the validated wide
accepted failure-start checkpoint (staged SHA-256
`0fd8e22e8ae04d51e829a1c2b19a0c8c71f8211f26de8a501f4c837da133ea46`)
and reached exact `t=2.55` in 20 accepted steps, native exit 0, wall 68.177 s.
Its [runner](kepler_sparse_clipcompact_cross_a10_20261027.slurm)
and staged inputs passed hash, shell-syntax, Python-import, and native-time
preflight. The postrun state/EOS, radial, native-versus-ballistic density,
and seam-provenance audits all passed. The separate full
[physical gate](KEPLER_SPARSE_WIDE_CLIPCOMPACT_PHYSICS_T255_2026-10-27.json)
also passed all 14 checks: finite positive state, periodic area, fixed mass,
W=0, `rho=m/V`, ideal-gas EOS, zero floor/dual/entropy injection, total
energy relative error `3.09937e-8`, and disk mean absolute log-entropy
`5.21595e-6`. The original disk density change is `1.37752e-4`, its
independent exact-orbit prediction is `1.37840e-4`, and their pointwise
mass-weighted residual is `6.83898e-7`. This proves the corrected native
operator crosses the **specific 48-corner abort** without visible early disk
degradation; it is not yet a t>6 result. The job writes to the new directory
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/to0p825/clipcompact_from2517696_a10`;
the previous data and binary are retained. The held original H200 duplicates
**437174/437175** were canceled because they contain the same failed CUDA
binary; no accepted state was deleted.

An independent [all-cell periodic Voronoi audit at the patched t=2.55
endpoint](KEPLER_SPARSE_WIDE_CLIPCOMPACT_T255_CORNERS_2026-10-27.json)
now strengthens the geometry check beyond the disk: every one of the
55,931 physical cells is bounded, the independent areas tile 36, the
maximum true corner count is **9**, no cell has 48 corners, and the maximum
native-versus-independent relative area difference is `5.84e-8`.

At `t=2.55`, 28 initially seam-near cells carrying fixed mass `0.00684143`
have centers inside the `r=2.35` guard, but none inside `r<2`. The center
gap to the original disk has shrunk to `0.210625`. Original `2.35<=r<4.25`
cells have mean absolute density change `0.4347`, pressure change `1.0483`,
and angular-momentum change `-0.3784`, while disk metrics remain close to
ballistic. This is a strong warning about the current periodic exterior,
which must be carried into every later interpretation.

The [saved local timestep audit](KEPLER_SPARSE_WIDE_DT_ORIGIN_T255_2026-10-27.json)
locates the smallest hydro-CFL candidate at initial seam-near cell ID
55285 (`r0=3.542`, current `r=2.236`), `dt_hydro=0.003125`; the six
smallest hydro candidates all start within `0.077` of a box face. At this
checkpoint the global accepted limit is still the innermost core orbital
bound `0.00193984`, not that hydro limit. Later transient accepted timestep
contractions therefore require another saved-state or in-step argmin audit
before a specific limiting cell can be named; the t=2.55 evidence already
shows the seam cells are the closest hydro constraint.

The passing `t=2.55` endpoint SHA-256 is
`503dad7e489fdbf86b8869d3297ddb343675d6d200e73423a59d95e8420fea46`.

An independent early-regression A10 replay **437280** used the same accepted
wide `t=.825` checkpoint and stopped at exact `t=.835` with the compacted
binary. It completed in 8.418 s and passed the state/EOS, radial,
exact-orbit geometry and seam auditors. Its complete native snapshot differs
in file SHA-256 from the original profiling run's `t=.835` dump, but an
ID-aligned comparison found **zero differing values** in every physical
particle field (`x`, `y`, `vx`, `vy`, mass, volume, internal energy, density,
pressure, and W). Thus the isolated clipping edit leaves this ordinary
early evolution state exactly unchanged; file-level metadata differs.

After checking the `t=2.55` physics, state, radial, geometric and seam receipts,
an independent continuation **437279** was staged and submitted on one A10
to exact `t=6.5` with the same fixed-mass W0 physics, RK4, halo 32, and
corrected CUDA binary. Its [runner](kepler_sparse_clipcompact_to6p5_a10_20261027.slurm)
checks every staged input SHA before evolving and then requires the full
physical gate, state/EOS, radial, native-versus-ballistic, and seam audits.
The independent run directory is
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/to0p825/clipcompact_to6p5_a10`.
As of submission there is **no t>6 endpoint**. Treat the outer periodic seam
as a separate physical-boundary attribution problem, not as a GPU-capacity
issue.

## Halo-32 failure at accepted t=2.755142729 (2026-10-27)

Job **437279 FAILED 1:0** after 10m23s at the next RK4 stage from accepted
step-start `t=2.7551427294639468`. The native log reports
`[GPU_STENCIL_EXPANSION] ... halo_clipped_cells=4` and then
`[GPU_TESS_INVALID] volume=0 velocity_fit=0 geometry=0 halo=4`. It saved
`khdiag.step001522.start_geometry.dat`; the failed stage itself was not
serialized. This is a distinct limit from the earlier 48-corner clipping
abort. The unchanged patched binary had no positive-state or injection
failure in the last accepted step.

The accepted-state [physical](KEPLER_SPARSE_WIDE_CLIPCOMPACT_PHYSICS_T275514_2026-10-27.json),
[state/EOS](KEPLER_SPARSE_WIDE_CLIPCOMPACT_T275514_STATE_2026-10-27.json),
[radial](KEPLER_SPARSE_WIDE_CLIPCOMPACT_T275514_RADIAL_2026-10-27.json),
and [native-versus-exact-orbit](KEPLER_SPARSE_WIDE_CLIPCOMPACT_T275514_BALLISTIC_2026-10-27.json)
audits pass. Disk mean absolute density change is `1.56915e-4` versus
exact-orbit `1.56997e-4`; the pointwise native-orbit residual is `7.77107e-7`.
Disk mean absolute log-entropy change is `5.85108e-6`, total relative energy
change `3.24093e-8`, fixed mass and W=0 exact, and all floor/dual/entropy
injections zero. These numbers establish a healthy **accepted state**, not
an accepted RK4 stage after it or survival to `t>6`.

An [independent periodic Voronoi halo reconstruction](KEPLER_SPARSE_WIDE_CLIPCOMPACT_T275514_HALO32_2026-10-27.json)
uses 164,049 physical/image sites. Its 2.5-unit image margin exceeds twice
the largest edge-cell circumradius by a factor 1.128. At the saved accepted
state, all 9,465 audited edge cells are still strictly inside the 32-cell
halo bound, but the worst four symmetric cells have ratio **0.999888** and
need `1.999753` units of the available `2.0`. Their independent areas agree
with native saved areas. A later RK stage crossing this bound is therefore
consistent with the native four-cell failure; the failed-stage geometry
itself remains unmeasured. The 36-cell width passes this **one saved-state**
bound, but is not a certified long-time width.

The [seam provenance](KEPLER_SPARSE_WIDE_CLIPCOMPACT_T275514_SEAM_2026-10-27.json)
shows 36 initially seam-near/outer cells, fixed mass `0.008764`, already
inside the `r=2.35` guard; none has a center inside the original `r<2` disk.
The outer-origin to disk-origin center gap is only `0.07882`.
Seam-tagged mean absolute density change is `5.734`, and its maximum inward
radial speed is `0.755`. The [local timestep audit](KEPLER_SPARSE_WIDE_CLIPCOMPACT_T275514_DT_ORIGIN_2026-10-27.json)
places the minimum hydro bound, `0.000979`, at initially outer cell ID 55324,
below the core orbital bound `0.001940`: the corrupt exterior now controls
global timestep cost. Wider padding alone would correct geometric
completeness while retaining the non-equilibrium periodic seam. The next
physical branch therefore needs an isolated/bounded-domain boundary
contract; see the [boundary audit](KEPLER_ISOLATED_BOUNDARY_GAP_2026-10-27.md).
The separate full-q2 H200 control remains in progress.
