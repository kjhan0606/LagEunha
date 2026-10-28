# Distributed pressure-seed response near t10

Built a zero-net-internal-energy seed on r<.6 from the actual geometric
compression signal -gamma Vdot/V. Subtracted its U-weighted mean and scaled
max absolute value to one. Four independent checkpoint copies used fractional
pressure changes +/-1e-5 and +/-5e-6. Pressure, total cell internal energy,
sound speed and entropy K were adjusted consistently by ratio, ratio, sqrt
ratio and ratio. The K offset was compiled from the original C type/layout
(624-byte record, stress.K at byte 600), rather than guessed. Position,
velocity and mass were unchanged. Net U perturbations are below 8.3e-22.

All four full-solver branches completed eight .0005 steps from
t10.0020667693911 to 10.006066769391104, finite with zero retries. The
same-binary unperturbed reference was reused. The positive quadratic
diagnostic sum m|delta v|^2 + sum V(delta P)^2/(gamma P), normalized by its
initial value, is 1.00122582. Kinetic and pressure terms are .001424405 and
.999801417. Half-amplitude response differs by 1.33e-7 in the corresponding
normalized norm. Even-in-seed pressure deviation is at the numerical floor.

Thus this distributed pressure/entropy seed persists over the tested short
interval while the earlier local velocity seed strongly damps. A .12% increase
in this norm is not established unstable growth: background weights evolve,
displacement/entropy terms are omitted, and no time-step convergence or
longer trajectory comparison was made for this seed. It is an initial
entropy perturbation, not a zero-entropy volume displacement, and is not the
same as continually forced geometric compression. No production steering or
pressure resetting occurred.

Next: continue preserved branch endpoints over a longer common interval and
compare with the prior velocity-seed response, testing step-size sensitivity
before interpreting a small excess above unity as physical amplification.

Evidence: canonical_pressure_seed_response_2026-10-18.json.
Scripts: run_pressure_seed_branches.py, analyze_pressure_seed_branches.py.
Remote outputs: /gpfs/kjhan/LagForce/Codex_pressure_seed_20261018.
