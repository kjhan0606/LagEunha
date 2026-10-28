# Weight-induced mass remap on a native Kepler step

The native `t=6.5→6.501` weighted experiment reduced fixed-mass cell-volume
change, but a power face that moves because its weights change sweeps fluid.
The [read-only midpoint audit](../audit_kepler_weight_only_mass_remap_20261027.py)
estimates that one geometric contribution to conservative mass transfer on
the saved aligned and 5%-jitter states. It uses the same hash-verified W0
and weighted starts, actual endpoint-average weight rates, independently
rebuilt power faces, arithmetic face density, and all faces touching the disk
(including disk/exterior faces). Each transfer is equal and opposite.

For a face from `i` to `j`, the weight-induced normal face speed is
`w_n^W=(Wdot_i-Wdot_j)/(2d)`. A density `rho_f` swept into the expanding
cell changes its mass by `rho_f L w_n^W dt`. The audit adds this hypothetical
mass change to the **saved weighted endpoint volume**; it does not evolve
the native solver.

| Initial disk mean absolute density change over the step | Aligned | 5% jitter |
|---|---:|---:|
| Native W0 fixed mass | `3.74859e-6` | `8.60583e-6` |
| Native weighted fixed mass | `4.29208e-7` | `4.10286e-7` |
| Weighted endpoint plus weight-only mass remap | `3.74601e-6` | `8.60329e-6` |

The remapped density-change field differs from W0 by mass-weighted mean
absolute `2.30e-8` aligned and `3.25e-8` jitter. The latter's p99 absolute
cell mass change is `5.40e-5` of its mass. Summed transfers cancel across
all touched cells to `1.7e-23` absolute or better; independent midpoint
power areas agree with midpoint native volumes within `2.95e-6` relative.
[Aligned receipt](aligned.json) · [Jitter receipt](jitter.json).

The result follows from the local continuity identity. If face density is
uniform, the weight-induced mass change is `dm_i^W=rho dV_i^W`, so its
contribution cancels in `d rho_i=(dm_i-rho_i dV_i)/V_i`. In this one-step
estimate, accounting for that swept mass almost removes the apparent density
benefit of holding mass fixed while changing power-cell volumes.

This is a midpoint, first-order geometric remap using endpoint-difference
weight rates. It is not an RK4 stage ledger, a full Riemann/ALE flux, or a
momentum/energy/spin-conservative evolution. It does not rule out every
dynamic-weight method or the separately queued fresh-from-`t=0` test. It
does show why a proposed conservative weight-linked mass exchange must be
tested against the density benefit instead of assuming both follow from the
fixed-mass pilot. The separate [W0 reduction analysis](../KEPLER_W0_ALE_FIXED_MASS_CONTRACT_2026-10-27.md)
records the moving-face consistency constraint.
