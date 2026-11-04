# Actual FV vortex/cold-flow qualification: cold gate failed

The native FV--ALE prototype published in `1a79410` improves the measured
Gresho vortex but fails the registered cold-flow gate. It is not adopted
as a common physical solver, and KH is not launched with this failed
candidate. This is evidence from actual native evolution, not rejection
based only on an initial projection estimate. The full user objective
remains active.

## Execution and physical reference

All five registered smooth cases use the unchanged v3 C/OpenMP binary,
SHA256 `4775b83a3533decf0791815c584f35da288b4a2386a3375bdff1a874c1514510`,
under `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/fv_native_smooth_v1`.
Two parallel slots request eight threads each on disjoint physical cores
0--7 and 8--15, verified from the current CPU topology. No state, limiter,
mu, quadrature, solver cap or physical operator was adjusted per test.

The new read-only `fv_smooth_cell_reference.c` integrates the stationary
analytic solution on each actual native power polygon, using positive
centroid triangles and order-8/order-16 Duffy quadrature controls.
`assess_fv_smooth.py` reports actual-cell mass-mean velocity and physical
mean pressure separately from reconstructed physical velocity at the
generator. The comparison neither subtracts a kinetic estimate from the
evolved state nor changes any checkpoint. Reference quadrature differences
are much smaller than the measured evolution errors.

## Gresho: genuine refinement improvement, not complete qualification

All three vortex runs reach their full t=1 endpoint without floors or
energy repairs. The initial pressure projection discrepancy from the
order-16 polygon reference is about 2.30e-9 in absolute volume L1 at 64^2.

| Run | Physical node velocity L1 | Cell-mean velocity L1 | Physical pressure relative L1 | Mean-state thermal energy increase |
|---|---:|---:|---:|---:|
| 64^2 | 0.0126674 | 0.0122853 | 0.000722157 | 0.000752864 |
| 64^2, step factor 1/2 | 0.0122797 | 0.0119160 | 0.000765453 | 0.000760839 |
| 128^2 | 0.00631061 | 0.00626891 | 0.000585069 | 0.000187772 |

The existing W0 baseline has physical node vector errors around 0.073 at
64^2 and 0.0381 at 128^2. The new physical node norms improve over these
values, and refinement nearly halves velocity error. The operators and
declared state projections differ; this is not an integrator-only
comparison. Initial analytic continuum fields and EOS are shared.

The cumulative absolute GCL volume defect at 64^2 decreases from
7.95857e-6 to 2.10633e-6 under step halving, a factor 3.78. This supports
temporal improvement of this diagnostic, not exact GCL. Vortex velocity
error changes by only about 3%, and thermal increase is nearly unchanged.
The observed spatial error and dissipation therefore remain material.

Evolution wall times are 54.26, 108.48 and 471.59 seconds, respectively.
These CPU diagnostics are not matched GIZMO timings. Geometry and weight
solve consume about 171.7 and 181.8 seconds of the 128^2 run, respectively;
4375 rejected attempts accompany 5644 accepted steps. Physics and cost
both require further work.

## Cold sine shear: endpoint alone gives a false impression

The exact solution has rho=1, ux=sin(2*pi*y), uy=0 and P=1e-6 for all
time. Initial relative pressure error is about 8--9e-12. The registered
criterion is full t=0.25, positivity and physical pressure relative L1
<=0.01, with refinement improvement.

| Run | Last accepted time | Endpoint reached | Pressure relative L1 at that time |
|---|---:|---|---:|
| 64^2 | 0.114772773 | no | 2.54558 |
| 128^2 | 0.25 | yes | 1.90519 |
| 128^2, step factor 1/2 | 0.127219979 | no | 0.671052 |
| 128^2, step factor 1/4 | 0.25 | yes | 1.29305 |

The additional time controls were registered as diagnostics after the
failure, not new adoption settings. Earlier failed endpoints are not
same-time error comparisons. No timestep setting is selected to hide a
failure. The two same-time 128^2 norms still differ substantially; pure
spatial-error dominance is not established. Nevertheless the quarter-step
result also misses 1% by more than two orders of magnitude.

The 64^2 run approaches the positive-pressure boundary, then encounters
inadmissible reconstructed fan states near floating-point cancellation.
The half-step 128^2 run also approaches zero pressure and fails. The
accepted pre-failure states and all diagnostics are retained; a retry or
solver-cap increase is not used to continue them.

Global total energy remains conserved to about 1e-16. The GCL energy
diagnostic divided by initial thermal energy is 1.86e-8 (64^2), 7.44e-7
(128^2), and 6.83e-8 (128^2 quarter-step). These diagnostics are small
beside the observed pressure error, but are not rigorous pressure bounds.

## Thermal chain rule separates conserved energy from heating

`fv_cold_rhs_audit.c` evaluates the unchanged native RHS on saved states:

    Udot = Edot - ubar dot pdot + 0.5*|ubar|^2*mdot.
    Pdot = (gamma-1)*(Udot-U*Vdot/V)/V.

At birth the native total thermal rates are +6.64135e-6 at 64^2 and
+4.16086e-7 at 128^2, approximately a factor 16 reduction. The total
energy rate is zero and Vdot is zero. An analytic physical face-flux
control on the same geometry and mesh velocity gives a thermal rate
consistent with zero. Thus initial projection noise or time integration
alone does not explain the first-RHS numerical heating.

At the accepted 128^2 t=0.25 state, native thermal rate is 6.58841e-6;
the analytic physical-flux proxy gives 6.41210e-6. Their pressure-rate
absolute L1 values divided by initial P are 14.5113 and 14.4835 per unit
time. This proxy applies the chain rule with the *current native mean
velocity*, rather than an independently integrated exact cell mean.
It suggests a substantial representation/transport contribution on
deformed cells, but is not an exact measurement of Kvar_dot and does not
prove impossibility of conservative FV methods.

Since Edot sums to zero, a conservative energy-flux-only replacement
cannot change the instantaneous global thermal rate when mass and
momentum rates are held fixed. The next closure must consider their
kinetic compatibility together. A floor, energy reset or hidden thermal
reservoir would conceal this problem rather than resolve it.

## Comparison-source check and restart scope

The retained pristine public GIZMO source at commit
`a828c4ba79d67093bd1a3d04f6f90b5d2d94d65f` has an MFM adiabatic pair-energy
face correction in `hydro/hydro_core_meshless.h` and conserved-to-thermal
chain-rule conversion in `hydro/hydro_toplevel.c`. The global
`ENERGY_ENTROPY_SWITCH_IS_ACTIVE` is commented out by default in
`allvars.h`; its presence in `kicks.c` is not evidence that it is enabled.
This source inspection is not a new matched GIZMO simulation or proof
that the MFM face prescription transfers to the present weight-linked
variable-mass power cells.

Separately, `fv_restart_replay_v1` checks accepted true-power 1D/2D/3D Sod
checkpoints. Loading and writing the initial state is byte-identical,
including birth mu/m0/W0, integrated residual and boundary ledgers. Two
single-threaded continuations from t=0.2 to t=0.201 have byte-identical
final checkpoints in each dimension. This qualifies basic same-binary
replay, not cross-platform, GPU or MPI restart compatibility.

## Evidence

Run and read-only analysis groups are preserved in
`Research/FVSmooth_20261105`; original executables and full data remain in
the authorized GPFS project run root. The same-time comparison figure
`fv_smooth_figures_v1/smooth_physics_gate.png` shows vortex refinement and
the cold pressure profile, including its large systematic error.
The grouped Astra evaluation selected the replacement cell-mean KEP core.
Its subsequent coupling and remaining stress error are documented in
`UNIFIED_HYDRO_KEP_NATIVE_2026-11-05.md`. Neither failed candidate is
promoted to the manuscript baseline.
