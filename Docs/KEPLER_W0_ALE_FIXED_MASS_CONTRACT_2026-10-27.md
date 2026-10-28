# W=0 fixed-mass reduction and exact moving-face flux

The approved reduction requires Voronoi geometry and fixed individual cell
mass when `W=0`. For a geometrically moving face with normal velocity `w_n`,
the exact Eulerian/ALE mass flux out of a cell is

`F_m = rho_f (v_f·n - w_n) L`.

For an incompressible uniform-density flow, exact mass inside the moving
polygon obeys `dm_i/dt = rho dV_i/dt`. A W0 Voronoi face is generally not
material even if every generator follows the fluid velocity. Consequently,
the exact geometric ALE flux need not preserve an individual Voronoi cell's
mass. The [small reproducible shear test](audit_kepler_w0_ale_fixed_mass_20261027.py)
uses `v=(0,x)`, constant `rho=1`, no viscosity and 49 unsteered generators.
For one bounded six-face cell the [receipt](KEPLER_W0_ALE_FIXED_MASS_SHEAR_2026-10-27.json)
records:

| Quantity | Value |
|---|---:|
| Cell area | `1.0646331285` |
| `dV/dt`, finite difference | `-5.208347542e-3` |
| `dV/dt`, sum of geometric face speeds | `-5.208347562e-3` |
| Sum of physical fluid normal flux | `1.60e-16` |
| Exact ALE `dm/dt` | `-5.208347562e-3` |

The finite-difference and face-speed rates agree to `2.0e-11`. Holding that
cell's mass fixed instead produces an instantaneous density change of
`+4.89215e-3` per unit time, although the exact uniform-density shear remains
uniform. This is a geometric counterexample, not a native Kepler error
estimate and not a claim that a fixed-mass discretization cannot converge.
The actual Kepler face-slip proxy likewise found substantial W0 face
crossings in the [one-step native checkpoint analysis](KEPLER_Q2_NATIVE_WEIGHTED_T6501_2026-10-27/README.md).

This establishes a design constraint: an operator that uses the **full
exact geometric ALE flux** for arbitrarily small nonzero weights cannot also
have a continuous `W→0` limit of exactly fixed individual masses for every
flow on a non-material Voronoi mesh. Total mass can remain pairwise
conserved in either case. The user's W0 fixed-mass reduction remains a
required discrete behavior; the extension must state explicitly how its
weight-linked mass exchange approaches that behavior and quantify the
remaining face-motion error under resolution refinement. A special W0
branch, a weight-linked correction to the finite-mass baseline, or a face
motion that is material in the relevant limit are distinct design choices;
none is silently assumed to be the full geometric ALE law.
The [native one-step weight-only remap estimate](KEPLER_WEIGHT_ONLY_MASS_REMAP_2026-10-27/README.md)
also shows that simply transferring the mass swept by weight motion nearly
removes the observed fixed-mass density benefit on that step.

For a future inviscid face-flux/spin branch, equal-and-opposite mass,
momentum and total-energy transfers must be recorded on each interior face.
If the momentum transfer is applied at face centroid `c`, its moment about
each generator supplies the corresponding spin exchange, with the generator
motion term included when generator velocity and momentum are not parallel.
Total energy must include bulk kinetic, rotational and thermal energy;
changes of cell inertia and mass must enter that same ledger. The
[current variational-pressure torque analysis](KEPLER_VARIATIONAL_SPIN_PAIR_TORQUE_2026-10-27/README.md)
shows why the present generator force alone cannot specify a unique spin
source. At `Re=∞`, this construction must not introduce a physical viscous
traction. A native conservative flux, positivity tests, W0 reduction and
long-time Kepler convergence are still required.
