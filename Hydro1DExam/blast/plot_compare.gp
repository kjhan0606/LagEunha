set terminal pngcairo size 1400,1200 enhanced font "Arial,13"
set output "blast/compare.png"

set multiplot layout 3,3 title "Interacting Blast Waves (N=200)" font ",16"

cRef   = "#999999"
cBase  = "#000000"
cAAKs  = "#e41a1c"
cSPH   = "#984ea3"

set style line 1 lc rgb cRef  lw 2 dt 1
set style line 2 lc rgb cBase lw 1.5 pt 7 ps 0.35
set style line 4 lc rgb cAAKs lw 1.5 pt 7 ps 0.35
set style line 5 lc rgb cSPH  lw 1.5 pt 7 ps 0.35

set xlabel "x"
set key top right font ",10" spacing 0.8

# Row 1: Density
set ylabel "ρ"

set title "Base (w=0)\nL_1=9.106e-2"
plot "blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "blast/base.dat" u 1:2 w lp ls 2 t "Base"

set title "AA+κ_s=15 (w=0.08)\nL_1=8.583e-2  (−5.7%)"
plot "blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "blast/aa_ks15.dat" u 1:2 w lp ls 4 t "AA+κ_s=15"

set title "SPH\nL_1=1.304e-1"
plot "blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "blast/sph.dat" u 1:2 w lp ls 5 t "SPH"

# Row 2: Velocity
set ylabel "v"

set title ""
plot "blast/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", \
     "blast/base.dat" u 1:3 w lp ls 2 t "Base"

set title ""
plot "blast/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", \
     "blast/aa_ks15.dat" u 1:3 w lp ls 4 t "AA+κ_s=15"

set title ""
plot "blast/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", \
     "blast/sph.dat" u 1:3 w lp ls 5 t "SPH"

# Row 3: Pressure
set ylabel "P"

set title ""
plot "blast/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", \
     "blast/base.dat" u 1:4 w lp ls 2 t "Base"

set title ""
plot "blast/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", \
     "blast/aa_ks15.dat" u 1:4 w lp ls 4 t "AA+κ_s=15"

set title ""
plot "blast/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", \
     "blast/sph.dat" u 1:4 w lp ls 5 t "SPH"

unset multiplot
