# Cell spin status and conservative gate for the native Kepler solver

## What the current H200 polar runs contain

The q1/q2/q2c 384² runs use the user's CUDA+MPI+OpenMP solver with unweighted
Voronoi cells (`W=0`), fixed cell masses after the one-time IC assignment,
RK4, and no physical shear viscosity (`Re parameter in Voro = 0` in the
native input). The active [variational W0 pressure operator](KEPLER_POLAR_OPERATOR_IDENTITY_2026-10-27.md)
uses geometric face work; computed HLLC star pressure and short-pair pressure
do not drive its force or work. The run scripts
unset `GFS_LAGUERRE_ROTATION`, `GFS_DYN_W2`, and `GFS_W2_LINEAR`; there is no
cell-spin state in the serialized particle layout or RK4 variables.

`GFS_LAGUERRE_ROTATION` is **geometric motion of a Laguerre face**, despite
its name. In `Exam/exam.c::voronoi_face_rotation` and the CUDA face kernel it
adjusts the face's normal speed for the face centroid's offset from its
weighted bisector anchor. The correction is already applied on W=0 faces
without enabling that flag. It stores no per-cell internal angular momentum,
moment of inertia, angular speed, or rotational energy. Enabling it in the
present W=0 run would therefore not test the user's cell-rotation proposal.

The earlier [native moment/spin audit](/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/RESULT_FAST_FV_COLD_MEAN_SPIN_20261026.md)
used a different 71,080-cell refined mesh. Its physically defined rigid
spin `S²/(2I)` represented 10.1677% of the internal velocity variance; the
central 52 cells still had a median residual variance 2.61639 times their
true thermal energy after spin was removed. This rules out scalar spin alone
as that layout's cold-EOS repair. Separately, the C010 W=0 disk's density
change at t=6.5 was almost entirely predicted by exact-orbit Voronoi area
change, which spin cannot alter when geometry and mass are fixed. The new
polar layout may remove that geometric error; native long-time results are
still pending.

## Minimum conservative cell-spin operator

For a shared face at centroid `c`, let the pressure/Riemann force on cell `i`
be `F` and on `j` be `-F`. If cell centers are `x_i,x_j`, orbital angular
momentum changes by `(x_i-x_j) × F dt`. Internal spin increments must be

`dS_i = (c-x_i) × F dt`, `dS_j = -(c-x_j) × F dt`.

Their sum cancels the orbital change pairwise. An off-center pressure force
can therefore transfer orbital and internal angular momentum without
inventing tangential friction. A Kepler run at physical Re=infinity must not
add a shear traction merely to damp numerical noise. If a finite Re is
specified later, the tangential stress and its positive heat production need
their own constitutive law and pairwise budget.

The center in this decomposition must be the **mass centroid**, with its
first-moment evolution included; it need not coincide with the tessellation
generator or the uniform-polygon area centroid. The newer
[centroid closure analysis](KEPLER_SPIN_CENTROID_CLOSURE_2026-10-27.md)
and [two-cell ledger](KEPLER_SPIN_FACE_LEDGER_PROTOTYPE_2026-10-27.md)
make this distinction explicit. Also, the variational generator force alone
does not uniquely specify a physical per-face traction or torque; the
face decomposition must be justified before using the above accounting as
a native spin law.

A useful native state would need internal spin `S_i`, polygon inertia `I_i`,
and rotational energy `S_i²/(2I_i)` with a consistent update when the cell
shape changes. Reconstructed velocity at the face must include the spin
contribution `omega_i z-hat × (c-x_i)`; otherwise the new state cannot affect
the Riemann shear error. Bulk kinetic, spin kinetic, thermal energy, and face
work must close one energy ledger. W=0/S=0 must recover the current Voronoi
solver. Polygon area, centroid, and second moment can in principle be
accumulated in O(number of faces), avoiding costly subcell quadrature; CUDA
and MPI halo exchange must carry the same moments and state.

## Decision gate

First compare the native q1/q2/q2c density, entropy, orbit and angular
momentum at the scheduled checkpoints. If polar alignment removes the disk
area noise but a residual shear/heating or angular-momentum error appears,
audit the face force's moment arm and the numerical tangential flux before
adding spin. A cell-spin branch needs a pairwise angular/energy audit, a cold
initial pressure gate, and a measured CUDA+MPI step cost before its Kepler
run is promoted beyond t=6. Spin is part of the full research objective,
but the current W=0 polar test does not yet contain it.
