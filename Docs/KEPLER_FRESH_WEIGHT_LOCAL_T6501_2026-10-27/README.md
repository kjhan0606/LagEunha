# Fresh-weighted Kepler: local work and swept mass at t=6.5

The accepted fresh-weighted `t=6.5` checkpoint was replayed for one RK4 step
to `t=6.501` on an idle RTX 5000 Ada on lageunha. It used the original v5
native executable, identical physical settings and accepted parent hash, in
an isolated directory. The [native endpoint receipt](audit_t6501.json) passes
all 21 checks; the [stage manifest](stage_manifest.json) records the input
hashes. The ongoing Slurm continuation was not modified.

CPU job **437575** reconstructed the 148,557 active-disk-touching power faces
at the midpoint using endpoint-secant weight rates. Its
[local pressure-work result](local_work.json) gives the 99th percentile of
`|P Vdot_weight| dt / U` as `3.65177e-5`; the total volume change has RMS
`2.09708e-7`, versus `1.41441e-5` for the weight-only component. The weight
and generator-motion terms mostly cancel in the total volume rate. The
[endpoint PdV check](pdv_closure.json) finds the 99th percentile of
`|ΔU + (P_start+P_end)ΔV/2|/U_start` is `8.46541e-8`. This is evidence of
smooth one-step cell thermodynamic consistency, not a four-stage pairwise
flux or material-surface proof.

CPU job **437578** applied a separate frozen-state, first-order geometric
estimate of mass swept by the weight-only face speed to the saved native
endpoint. This is not an evolved mass flux. The disk mass-weighted mean
absolute one-step density changes are:

| Estimate | Mean absolute relative density change |
|---|---:|
| Original fixed-mass weighted step | `1.28250e-7` |
| Add sweep using arithmetic face density | `8.68398e-6` |
| Add sweep using weight-sweep upwind density | `8.68399e-6` |

Both mass estimates are pairwise conservative and positive for this step;
their p99 difference in individual mass change is `9.95e-8` of cell mass.
Midpoint reconstructed power area agrees with averaged native area within
`2.92e-6` relative. [Result](weight_sweep.json).

The experiment is on a fresh weighted trajectory, rather than imposing
weights only after an unweighted history. It confirms that the geometry-only
mass exchange would substantially change the observed one-step density
response. It does not compare against a same-start W0 step, prescribe the
full ALE upwind flux, implement momentum/energy transport, or establish the
long-term result of a self-consistent variable-mass evolution. A native RK4
flux and physical validation are still needed before a conservation claim.
