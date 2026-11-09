# Common bounded-velocity comparison, 2026-11-09

## Historical reconciliation and decision

The preceding common-mode proposal overlooked a documented rejection. The
unbounded common-speed projection was already tested in `D/continuous_jump`
and repeated in `S/continuous_normal_v2`; see
`RESEARCH_PRIORITY_AND_EVIDENCE_2026-11-08.md`. It must not be repeated again
as a new idea. The recent `D/common_velocity_pair_20261109_v1` local tests
and 1D comparison are preserved but rejected. With the current scalar-v3
reconstruction its nine 1D cases reach their endpoints, yet Blast density L1
is 0.09929531 versus scalar-only 0.04680570 and original 0.04110979;
colliding shocks worsen to 0.00803596. No native 2D evolution or adoption of
that unbounded candidate follows. Completion and algebraic conservation do
not override this evidence.

Here D is `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1`.
The active isolated experiment is `D/common_bounded_velocity_20261109_v1`.
It reuses `D/bounded_jump/gfs_material_bounded_jump.h` byte-for-byte
(SHA256 34a99d6f5bbb2690800932643225b083dc3d33959fa27bd2aa9218e60d8cec76).
This historical helper projects the two normal traces jointly onto their
stencil bounds and admissible relative-jump strip. Its mathematical formula,
physical variables and coefficients are unchanged. The new combination is
the scalar-v3 reconstruction and the current direct-FV C/CUDA adapter, not
a newly invented velocity limiter. It is not adopted as the default.

## Implementation and necessity

The previous alpha multiplication can suppress a finite common face speed
when the owner normal-speed difference is zero. The stored irregular-cell
affine fixture reproduces that failure. The bounded helper repairs that
specific example without allowing the unlimited common speed that fails
shocks. Raw traces, force/volume/energy exchange, and cached pressure-response
states use the same paired reconstruction. No new stencil entries, physical
state, empirical coefficient or test-specific sensor is introduced.

The adapter initially computed stencil projections with a two-term dot
expression, whereas the historical helper accumulated components from zero.
One CUDA fixture rejected its own owner bound because of arithmetic ordering.
Using the identical accumulation in the adapter resolves that fixture without
expanding bounds or adding a tolerance. Original failure and source recipe
are retained. `recipe_reconciliation.json` verifies that the updated recipe
emits the already tested adapter exactly. The historical helper is unchanged.

This candidate is restricted to the periodic direct-FV path. Prescribed
exterior/wall traces are explicitly unsupported; no cylinder/RT/Kepler result
is attributed to it. Frozen production binaries and running states are intact.

## Completed checks and their scope

- Actual local affine face-speed error is 1.73e-17; perturbations of either
  sign at 1e-14 and 1e-10 produce correspondingly small trace changes.
- The recorded vanishing-face wave fixtures have force changes decreasing
  approximately with geometric perturbation, down to about 7e-9 at 1e-8.
  This is fixture evidence, not a proof for arbitrary meshes.
- Eleven complete C/CUDA RHS fixtures agree near floating-point precision.
  Nine wave fixtures also test the exact CUDA object linked into the native
  executable, with maximum scaled discrepancy 1.83e-15. The latter avoids
  treating a separately optimized standalone test as runtime proof.
- All nine 1D cases finish, with fixed mass and energy-balance residuals near
  roundoff. Density L1: Sod .00402172, Blast .05317852, Noh .00408063,
  Shu .00296239, Lax .00301232, double rarefaction .02312789, colliding shocks
  .00398906; contacts remain at roundoff. Blast remains worse than both
  original and scalar-only controls. These are mixed accuracy results.
- Eighteen actual initial/final 1D states embedded in weighted rectangular
  strips agree with the 2D semi-discrete interior operator within 4.6e-14;
  transverse force is zero. This excludes differing outer boundary stencils
  and does not certify full moving-geometry evolution, MPI or 3D.

Evidence: `adapter_continuity.json`, `cuda_report.json`,
`native/runtime_object_parity.json`, `one_d/candidate_accuracy.json`,
`dimensional_reduction/report.json` and their retained sources/inputs.

## Native comparisons underway

Native executable SHA256:
`c5e25064a2bada3fad8b3119f645f83f6f89ae3887719ebd04f99e598139a2ad`.
Only exam_gpu_extract.o and material_cuda_rhs.o differ from the scalar-only
parent archive. Both 32-square axis/diagonal wave endpoints and actual-cell
assessments have completed. The predeclared sequential 64-square shear and
Gresho comparisons then began, with the same archived ICs and settings.
Controls are reused; they are not rerun. All segments continue their own
checkpoints. See `native_batch_plan.json`, `waves/status.json`, and
`flows/results.json` for authoritative endpoint state. A live process was
verified, rather than relying on a status file to infer execution.

The current combination is not adopted. Local continuity, endpoint completion
or reduced heat cannot replace full error and shape comparisons. Boundaries,
full second-order convergence, MPI, 3D, and end-to-end <=2x GIZMO at matched
accuracy remain unproven. Finish 2D before starting 3D. The first astronomy
draft retains priority; no second manuscript is begun by these diagnostics.

## Completed acoustic/shear comparison

The independent read-only audit and `completed_wave_shear/report.json` confirm
byte-identical IC checkpoints, extracted IC fields and parameter files for
the three acoustic methods. Same actual-cell quadrature gives:

|32-square, one period|Original density L1/epsilon|Scalar only|Scalar + bounded velocity|
|---|---:|---:|---:|
|Axis|.04003255868|.01951583696|.01975259836|
|Diagonal|.05702601917|.01991670838|.02002552513|

The bounded velocity addition raises density error by1.213%/0.546% relative
to scalar-only. Most improvement against original control was already due to
the scalar change. These are two directions at one resolution, not a
convergence study. The stale 'Four cases' scope wording in the bounded wave
plan is corrected to 'Two cases'; the original manifest remains preserved.

All four shear segments completed to1. Thermal-energy increase is
0.000484235046016 relative to initial U (0.0484235%), essentially identical
to scalar-only. The maximum absolute difference across the13 extracted
state fields is5.063e-16. Thus the current smooth-shear test supplies no
evidence of additional accuracy gain from the bounded velocity change.
Masses remain fixed, U/V positive, energy relative error below4.6e-15.

Gresho remains in progress toward its predeclared1.2 endpoint. A finite
postprocessor verifies the live supervisor identity while waiting, then
compares frozen original, scalar-only and bounded methods with the unchanged
actual-cell quadrature and thermal budgets. It writes no checkpoints and
never restarts a simulation. See `postprocess/plan.json` and status.
The final shape report must be inspected before any method decision.

The cached direct pressure-response state was inspected in source and uses the same bounded traces and projected thermodynamic values. An independently exercised active-constraint trajectory has NOT been certified by these periodic runs; source consistency alone is not that runtime gate.

## Final Gresho result and strong-shock rejection

All nine flow segments and the saved-state postprocessor completed. Actual
polygon images were inspected with common scales and no smoothing. At1.2:

|Method|Acoustic heat|Cell momentum L1|Cell pressure L1|Density L1|
|---|---:|---:|---:|---:|
|Original|.003813917241|.04082231370|.001013465140|.001034765641|
|Scalar only|.003806064800|.04315202791|.001039033937|.001050320556|
|Scalar + bounded velocity|.003225071907|.04019206409|.000817875902|.000944793936|

Thus this final vortex test does show15.44% less acoustic heat and19.30% less
pressure error than original; momentum error improves only1.54%. The earlier
no-additional-benefit statement applies to wave/shear, not this endpoint.
The scalar-only and unbounded candidates must not be conflated with this one.

A predeclared two-layout64-square Sedov batch then used this SAME executable,
coefficients, original IC checkpoints, and original parameters. Glass seed42
was first. It exited99 after one accepted step at2.6856168702688705e-6; the
latest ordinary saved checkpoint is still0. The second regular-grid case was
NOT started, following the declared failure stop. No trial was restarted or
retuned; old successful glass and regular controls remain intact.

The unchanged C face-response replay reproduces native CPU/GPU status2:
reconstructed face pressure9.094947005735e-19, normal separation4.367136e-6,
vacuum threshold7.387148e-9, ratio591.18. Status2 means a vacuum fan has no
single shared material contact. This is not a CUDA-only discrepancy. The
pressure numerically matches the scalar margin4096epsilon times approximately
1e-6, but original owner/raw trace attribution is not yet reconstructed, so
that link remains a hypothesis. The diagnostic stop state is not promoted to
an accepted restart. See `sedov/glass/vacuum_audit/report.json` and `replay.c`.

**Decision: do not adopt this configuration.** Local continuity and better
Gresho heat do not suffice when the same candidate fails a strong2D shock.
Next investigate coupled thermodynamic/velocity trace admissibility using
this retained counterexample and prior positive-fan history. No added pressure
floor, test-specific coefficient, tolerance change or new3D campaign follows.


## Recovered thermodynamic origin of the Sedov counterexample

The prior paragraph's owner/raw-trace uncertainty is now resolved by
`sedov/glass/state_trace_audit`. Independent exact periodic power geometry
reproduces stored volumes within2.40e-14 and matches reverse faces within8.96e-16.
The production reconstruction is reused with diagnostic printing only.

Both owner pressures are healthy: 9.9999999868e-7 and1.0016258789e-6.
The right raw LS face pressure is -1.8750563 and its tangentially shifted value
is -0.138119741; the left raw face pressure is1.00026437e-6. The local cold owner
has a strong hot-neighbor gradient, rather than being an actual vacuum cell.
The scalar-v3 interval then projects both face pressures to9.0949470057e-19,
its4096epsilon margin scale. The resulting tiny acoustic speed produces the
already reproduced vacuum status2. No floor/tolerance or stored state was changed.
The stop diagnostic remains distinct from an accepted restart checkpoint.

A subsequent actual-stencil primitive interval removes this one static failure,
but a stronger disappearing-neighbor test falsifies that alternative: it leaves
a finite22.8% force jump as geometric perturbations tend to zero. See
`COMMON_STENCIL_PRIMITIVE_2026-11-09.md` and the candidate's `decision.json`.
Neither configuration is adopted; passing isolated states is insufficient.


## Additional velocity topology counterexample

The stronger fixed-velocity test under
`common_log_acceleration_probe_20261109_v1/continuity_velocity` also evaluates
this unchanged bounded-velocity parent. With rho=P=1, fixed per-cell velocities
in[-.1,.1], and vanishing geometric shear, its force jump approaches8.10% rather
than zero. The log-scalar candidate gives identical output in this constant
thermodynamic case. This independently rejects the inherited actual-stencil
normal-speed bounds as a generally continuous operator. The earlier smooth
wave fixture pass remains valid only for its narrower inputs. See
`COMMON_LOG_PRIMITIVE_2026-11-09.md` for the complete retained evidence.
