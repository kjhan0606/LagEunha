# KH native weight comparison through t=0.3

Existing native C/MPI/CUDA executable, Syntax A10, one MPI process and two
OpenMP threads, 16x32 cells in periodic 1x2 domain. Same smooth KH generator
and native adiabatic pressure force; Re infinite, fixed cell masses, no
steering or shocks. Dynamic resumes the preserved t=.1 snapshot; frozen is
its own newly required control trajectory from t=0. No dynamic replay.
Artifact: Research/MaterialVolume_20261105/kh_weight_comparison_20261106_v1.

Both reach .3, all subprocesses terminal rc0. Restart intervals <=.015
preserve the existing .02 window gate; each source checksum/environment is
recorded. Native timestep bounds differ between frozen and dynamic, and the
restart endpoints differ; this is not an identical-timestep comparison.

| Quantity at .3 | Frozen W0 | Dynamic W |
|---|---:|---:|
| max cell mass change | 0 | 0 |
| relative total energy error | 6.15e-10 | 7.21e-11 |
| max entropy invariant change | 2.39e-15 | 5.73e-15 |
| total area | 2 | 2 |
| pressure min/max | 9.97985 / 10.02029 | 9.96994 / 10.03105 |
| transverse velocity RMS | .00369615 | .00348031 |
| projected initial-mode amplitude | .00755148 | .00712696 |
| weight kinetic / translational energy | 0 | 3.91e-6 |

Projected amplitude is a volume-weighted least-squares coefficient of the
initial sin(2pi x) times two Gaussian shear-layer envelope, evaluated at
current generator positions. It is not a calibrated KH eigenmode growth
rate. Initial amplitude .01 need not grow monotonically for this initial
perturbation. No high-resolution/reference solution is provided here, so
neither candidate is declared more physically accurate from these numbers.
The coarse transition layer and short interval do not qualify nonlinear KH.

Dynamic weights again have no demonstrated universal advantage. Most .015
segments require ~39 dynamic steps versus 3--5 frozen steps; process startup
and repeated small windows dominate wall time here. Do not report this as
an end-to-end GIZMO benchmark. Next required work is longer KH with spatial
convergence, RT, common shock/entropy behavior with the successful 1D method,
and realistic runtime measurement. 3D remains deferred until 2D completion.
