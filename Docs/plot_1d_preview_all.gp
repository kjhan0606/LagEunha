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

# Velocity preview
set terminal pngcairo size 1400,1400 enhanced font "Arial,12"
set output "1d_preview_vel.png"
set multiplot layout 4,4 margins 0.08,0.97,0.06,0.95 spacing 0.055,0.04

set xrange [0:1]; set yrange [-0.05:1.05]
set ylabel "Sod  v" font ",13"; set title "SPH" font ",13"
plot "../sod/ref.dat" u 1:3 w l ls 1 t "Exact", "../sod/sph.dat" u 1:3 w p ls 2 t "SPH"
set ylabel "v" font ",13"; set title "Voronoi" font ",13"
plot "../sod/ref.dat" u 1:3 w l ls 1 t "Exact", "../sod/voronoi.dat" u 1:3 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot "../sod/ref.dat" u 1:3 w l ls 1 t "Exact", "../sod/price.dat" u 1:3 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot "../sod/ref.dat" u 1:3 w l ls 1 t "Exact", "../sod/hllc_cd10.dat" u 1:3 w p ls 5 t "HLLC+CD10"

set title ""; set xrange [0:1]; set yrange [-1:16]
set ylabel "Blast  v" font ",13"
plot "../blast/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", "../blast/sph.dat" u 1:3 w p ls 2 t "SPH"
set ylabel "v" font ",13"
plot "../blast/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", "../blast/voronoi.dat" u 1:3 w p ls 3 t "Voronoi"
plot "../blast/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", "../blast/price.dat" u 1:3 w p ls 4 t "Price AV"
plot "../blast/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", "../blast/hllc_cd10.dat" u 1:3 w p ls 5 t "HLLC+CD10"

set xrange [0:0.45]; set yrange [-1.1:0.15]
set ylabel "Noh  v" font ",13"
plot "../noh/ref.dat" u 1:3 w l ls 1 t "Exact", "../noh/sph.dat" u 1:3 w p ls 2 t "SPH"
set ylabel "v" font ",13"
plot "../noh/ref.dat" u 1:3 w l ls 1 t "Exact", "../noh/voronoi.dat" u 1:3 w p ls 3 t "Voronoi"
plot "../noh/ref.dat" u 1:3 w l ls 1 t "Exact", "../noh/price.dat" u 1:3 w p ls 4 t "Price AV"
plot "../noh/ref.dat" u 1:3 w l ls 1 t "Exact", "../noh/hllc_cd10.dat" u 1:3 w p ls 5 t "HLLC+CD10"

set xrange [-5:5]; set yrange [-4.5:3.5]
set xlabel "x" font ",12"; set ylabel "Shu-Osher  v" font ",13"
plot "../shuosher/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", "../shuosher/sph.dat" u 1:3 w p ls 2 t "SPH"
set ylabel "v" font ",13"
plot "../shuosher/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", "../shuosher/voronoi.dat" u 1:3 w p ls 3 t "Voronoi"
plot "../shuosher/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", "../shuosher/price.dat" u 1:3 w p ls 4 t "Price AV"
plot "../shuosher/ref.dat" u 1:3 w l ls 1 t "Ref(4000)", "../shuosher/hllc_cd10.dat" u 1:3 w p ls 5 t "HLLC+CD10"
unset multiplot

# Pressure preview
set output "1d_preview_prs.png"
set xlabel ""
set multiplot layout 4,4 margins 0.08,0.97,0.06,0.95 spacing 0.055,0.04

set xrange [0:1]; set yrange [0.05:1.1]
set ylabel "Sod  P" font ",13"; set title "SPH" font ",13"
plot "../sod/ref.dat" u 1:4 w l ls 1 t "Exact", "../sod/sph.dat" u 1:4 w p ls 2 t "SPH"
set ylabel "P" font ",13"; set title "Voronoi" font ",13"
plot "../sod/ref.dat" u 1:4 w l ls 1 t "Exact", "../sod/voronoi.dat" u 1:4 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot "../sod/ref.dat" u 1:4 w l ls 1 t "Exact", "../sod/price.dat" u 1:4 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot "../sod/ref.dat" u 1:4 w l ls 1 t "Exact", "../sod/hllc_cd10.dat" u 1:4 w p ls 5 t "HLLC+CD10"

set title ""; set xrange [0:1]; set yrange [0:460]
set ylabel "Blast  P" font ",13"
plot "../blast/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", "../blast/sph.dat" u 1:4 w p ls 2 t "SPH"
set ylabel "P" font ",13"
plot "../blast/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", "../blast/voronoi.dat" u 1:4 w p ls 3 t "Voronoi"
plot "../blast/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", "../blast/price.dat" u 1:4 w p ls 4 t "Price AV"
plot "../blast/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", "../blast/hllc_cd10.dat" u 1:4 w p ls 5 t "HLLC+CD10"

set xrange [0:0.45]; set yrange [-0.05:2.0]
set ylabel "Noh  P" font ",13"
plot "../noh/ref.dat" u 1:4 w l ls 1 t "Exact", "../noh/sph.dat" u 1:4 w p ls 2 t "SPH"
set ylabel "P" font ",13"
plot "../noh/ref.dat" u 1:4 w l ls 1 t "Exact", "../noh/voronoi.dat" u 1:4 w p ls 3 t "Voronoi"
plot "../noh/ref.dat" u 1:4 w l ls 1 t "Exact", "../noh/price.dat" u 1:4 w p ls 4 t "Price AV"
plot "../noh/ref.dat" u 1:4 w l ls 1 t "Exact", "../noh/hllc_cd10.dat" u 1:4 w p ls 5 t "HLLC+CD10"

set xrange [-5:5]; set yrange [0:10]
set xlabel "x" font ",12"; set ylabel "Shu-Osher  P" font ",13"
plot "../shuosher/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", "../shuosher/sph.dat" u 1:4 w p ls 2 t "SPH"
set ylabel "P" font ",13"
plot "../shuosher/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", "../shuosher/voronoi.dat" u 1:4 w p ls 3 t "Voronoi"
plot "../shuosher/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", "../shuosher/price.dat" u 1:4 w p ls 4 t "Price AV"
plot "../shuosher/ref.dat" u 1:4 w l ls 1 t "Ref(4000)", "../shuosher/hllc_cd10.dat" u 1:4 w p ls 5 t "HLLC+CD10"
unset multiplot

# Internal energy preview
set output "1d_preview_eng.png"
set xlabel ""
set multiplot layout 4,4 margins 0.08,0.97,0.06,0.95 spacing 0.055,0.04

set xrange [0:1]; set yrange [1.5:3.2]
set ylabel "Sod  e" font ",13"; set title "SPH" font ",13"
plot "../sod/ref.dat" u 1:5 w l ls 1 t "Exact", "../sod/sph.dat" u 1:5 w p ls 2 t "SPH"
set ylabel "e" font ",13"; set title "Voronoi" font ",13"
plot "../sod/ref.dat" u 1:5 w l ls 1 t "Exact", "../sod/voronoi.dat" u 1:5 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot "../sod/ref.dat" u 1:5 w l ls 1 t "Exact", "../sod/price.dat" u 1:5 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot "../sod/ref.dat" u 1:5 w l ls 1 t "Exact", "../sod/hllc_cd10.dat" u 1:5 w p ls 5 t "HLLC+CD10"

set title ""; set xrange [0:1]; set yrange [0:1500]
set ylabel "Blast  e" font ",13"
plot "../blast/ref.dat" u 1:5 w l ls 1 t "Ref(4000)", "../blast/sph.dat" u 1:5 w p ls 2 t "SPH"
set ylabel "e" font ",13"
plot "../blast/ref.dat" u 1:5 w l ls 1 t "Ref(4000)", "../blast/voronoi.dat" u 1:5 w p ls 3 t "Voronoi"
plot "../blast/ref.dat" u 1:5 w l ls 1 t "Ref(4000)", "../blast/price.dat" u 1:5 w p ls 4 t "Price AV"
plot "../blast/ref.dat" u 1:5 w l ls 1 t "Ref(4000)", "../blast/hllc_cd10.dat" u 1:5 w p ls 5 t "HLLC+CD10"

set xrange [0:0.45]; set yrange [-0.02:1.1]
set ylabel "Noh  e" font ",13"
plot "../noh/ref.dat" u 1:5 w l ls 1 t "Exact", "../noh/sph.dat" u 1:5 w p ls 2 t "SPH"
set ylabel "e" font ",13"
plot "../noh/ref.dat" u 1:5 w l ls 1 t "Exact", "../noh/voronoi.dat" u 1:5 w p ls 3 t "Voronoi"
plot "../noh/ref.dat" u 1:5 w l ls 1 t "Exact", "../noh/price.dat" u 1:5 w p ls 4 t "Price AV"
plot "../noh/ref.dat" u 1:5 w l ls 1 t "Exact", "../noh/hllc_cd10.dat" u 1:5 w p ls 5 t "HLLC+CD10"

set xrange [-5:5]; set yrange [0:13]
set xlabel "x" font ",12"; set ylabel "Shu-Osher  e" font ",13"
plot "../shuosher/ref.dat" u 1:5 w l ls 1 t "Ref(4000)", "../shuosher/sph.dat" u 1:5 w p ls 2 t "SPH"
set ylabel "e" font ",13"
plot "../shuosher/ref.dat" u 1:5 w l ls 1 t "Ref(4000)", "../shuosher/voronoi.dat" u 1:5 w p ls 3 t "Voronoi"
plot "../shuosher/ref.dat" u 1:5 w l ls 1 t "Ref(4000)", "../shuosher/price.dat" u 1:5 w p ls 4 t "Price AV"
plot "../shuosher/ref.dat" u 1:5 w l ls 1 t "Ref(4000)", "../shuosher/hllc_cd10.dat" u 1:5 w p ls 5 t "HLLC+CD10"
unset multiplot
