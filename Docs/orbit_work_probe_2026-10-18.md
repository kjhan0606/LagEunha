# RK circular-departure work instrumentation

Added opt-in GFS_ORBIT_WORK_AUDIT to Exam/exam.c. Captures the actual
TCF normal dissipative force per persistent ID, then reports each RK
stage's global active-particle hydro-minus-TCF, TCF and differential-rotation
work. For u=v-vc(x), the latter is -sum m*r*Omega'(r)*u_r*u_phi.
Uses the configured spline gravity analytically. Restricted to MPI1,
Kepler spline gravity, split HG, hydro enabled and zero XSPH/background
acceleration/Hill terms. Default-off does not change evolution formulas.

Logs are trial-stage records. Parsers must exclude rejected attempts;
current probe rejects any run with retries. Retry and MPI support for
this diagnostic have not been validated. Hydro-minus-TCF is not labelled
pure pressure when other hydro force options might be active.

Built a separate binary at Codex_orbit_work_20261018. Initial compilation
failed for a missing forward declaration of kepler_center; fixed before
successful build. Original running binaries and checkpoints are untouched.

From preserved t=2.1029059177634006 checkpoint, on/off ran one dt=.0025
step; half ran two dt=.00125 steps to the same endpoint within rounding.
All exited0, finite, no retries. On/off reader fields id,x,y,vx,vy,mass,
volume,den,pressure,ie,w2 are bitwise equal; unlisted fields were not checked.

| Increment over .0025 | Full step | Two half steps |
|---|---:|---:|
| Hydro minus TCF | +3.4722973e-12 | +3.4722722e-12 |
| TCF | -1.9664037e-12 | -1.9663669e-12 |
| Differential rotation | +1.5939777e-12 | +1.5951153e-12 |
| HG, directly measured | -8.2678281e-13 | -8.2675988e-13 |
| Actual departure-energy change | +2.2742559e-12 | +2.2742974e-12 |
| Endpoint minus HG minus RK quadrature | +1.1673976e-15 | +3.6605758e-17 |

The closure residual decreases by about32; this is one interval and not
a global formal-order claim. TCF and HG directly damp departure here;
hydro-minus-TCF and differential rotation supply it. This does not exclude
indirect dissipation-heating feedback. It does not establish long-time
stability, nor can it replace the historical t6-10 physical failure gate.

Next: extend the accepted-stage budget over a resolved early-growth interval
from the saved endpoint, retain closure residuals, and verify canonical
hydro-minus-TCF contains only intended pressure terms before naming it so.

Evidence: canonical_orbit_work_probe_2026-10-18.json.
Scripts: run_orbit_work_probe.py, analyze_orbit_work_probe.py.
