set terminal pdfcairo size 10,10 enhanced font "Arial,12"
set output "1d_tests.pdf"

set multiplot layout 4,4 margins 0.08,0.97,0.06,0.95 spacing 0.055,0.04

cRef   = "#aaaaaa"
cSPH   = "#984ea3"
cVoro  = "#4daf4a"
cPrice = "#377eb8"
cHLLC  = "#e41a1c"

set style line 1 lc rgb cRef   lw 3.0 dt 1
set style line 2 lc rgb cSPH   lw 0.5 pt 2 ps 0.7
set style line 3 lc rgb cVoro  lw 0.5 pt 2 ps 0.7
set style line 4 lc rgb cPrice lw 0.5 pt 2 ps 0.7
set style line 5 lc rgb cHLLC  lw 0.5 pt 2 ps 0.7

set key top right font ",10" spacing 0.8

# ============================================================
#  Row 1: Sod  (density)   rho: 0.11 -- 1.01
# ============================================================
set xrange [0:1]
set yrange [0.05:1.1]

set ylabel "Sod  {/Symbol r}" font ",13"
set title "SPH" font ",13"
plot "../sod/ref.dat" u 1:2 w l ls 1 t "Exact", \
     "../sod/sph.dat" u 1:2 w p ls 2 t "SPH"

set ylabel "{/Symbol r}" font ",13"
set title "Voronoi" font ",13"
plot "../sod/ref.dat" u 1:2 w l ls 1 t "Exact", \
     "../sod/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"

set title "Price AV" font ",13"
plot "../sod/ref.dat" u 1:2 w l ls 1 t "Exact", \
     "../sod/price.dat" u 1:2 w p ls 4 t "Price AV"

set title "HLLC+CD10" font ",13"
plot "../sod/ref.dat" u 1:2 w l ls 1 t "Exact", \
     "../sod/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"

# ============================================================
#  Row 2: Blast  (density)  rho: 0 -- 13
# ============================================================
set title ""
set xrange [0:1]
set yrange [0:13]

set ylabel "Blast  {/Symbol r}" font ",13"
plot "../blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "../blast/sph.dat" u 1:2 w p ls 2 t "SPH"

set ylabel "{/Symbol r}" font ",13"
plot "../blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "../blast/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"

plot "../blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "../blast/price.dat" u 1:2 w p ls 4 t "Price AV"

plot "../blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "../blast/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"

# ============================================================
#  Row 3: Noh  (density)   rho: 0 -- 10
# ============================================================
set xrange [0:0.45]
set yrange [0:10]

set ylabel "Noh  {/Symbol r}" font ",13"
plot "../noh/ref.dat" u 1:2 w l ls 1 t "Exact", \
     "../noh/sph.dat" u 1:2 w p ls 2 t "SPH"

set ylabel "{/Symbol r}" font ",13"
plot "../noh/ref.dat" u 1:2 w l ls 1 t "Exact", \
     "../noh/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"

plot "../noh/ref.dat" u 1:2 w l ls 1 t "Exact", \
     "../noh/price.dat" u 1:2 w p ls 4 t "Price AV"

plot "../noh/ref.dat" u 1:2 w l ls 1 t "Exact", \
     "../noh/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"

# ============================================================
#  Row 4: Shu-Osher  (density)  rho: 0 -- 4.5
# ============================================================
set xrange [-5:5]
set yrange [0:4.8]

set xlabel "x" font ",12"
set ylabel "Shu-Osher  {/Symbol r}" font ",13"
plot "../shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "../shuosher/sph.dat" u 1:2 w p ls 2 t "SPH"

set ylabel "{/Symbol r}" font ",13"
plot "../shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "../shuosher/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"

plot "../shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "../shuosher/price.dat" u 1:2 w p ls 4 t "Price AV"

plot "../shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", \
     "../shuosher/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"

unset multiplot
