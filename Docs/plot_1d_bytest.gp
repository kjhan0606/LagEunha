# ============================================================
#  Per-test figures: each PDF = 4 rows (rho,v,P,e) x 5 cols (methods)
#  Methods: SPH, Voronoi, Price AV, HLLC+CD10, NS+AV
# ============================================================

cRef   = "#aaaaaa"
cSPH   = "#984ea3"
cVoro  = "#4daf4a"
cPrice = "#377eb8"
cHLLC  = "#e41a1c"
cNSAV  = "#ff7f00"

set style line 1 lc rgb cRef   lw 3.0 dt 1
set style line 2 lc rgb cSPH   lw 0.5 pt 2 ps 0.7
set style line 3 lc rgb cVoro  lw 0.5 pt 2 ps 0.7
set style line 4 lc rgb cPrice lw 0.5 pt 2 ps 0.7
set style line 5 lc rgb cHLLC  lw 0.5 pt 2 ps 0.7
set style line 6 lc rgb cNSAV  lw 0.5 pt 2 ps 0.7

set key top right font ",10" spacing 0.8

# data directory relative to Docs/
D = "../src/Hydro1DExam"

# ============================================================
#  Sod Shock Tube
# ============================================================
set terminal pdfcairo size 12.5,10 enhanced font "Arial,12"
set output "1d_sod.pdf"
set multiplot layout 4,5 margins 0.07,0.97,0.06,0.95 spacing 0.045,0.04

# Row 1: Density
set xrange [0:1]; set yrange [0.05:1.1]
set ylabel "{/Symbol r}" font ",13"
set title "SPH" font ",13"
plot D."/sod/ref.dat" u 1:2 w l ls 1 t "Exact", D."/sod/sph.dat" u 1:2 w p ls 2 t "SPH"
set title "Voronoi" font ",13"
plot D."/sod/ref.dat" u 1:2 w l ls 1 t "Exact", D."/sod/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot D."/sod/ref.dat" u 1:2 w l ls 1 t "Exact", D."/sod/price.dat" u 1:2 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot D."/sod/ref.dat" u 1:2 w l ls 1 t "Exact", D."/sod/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"
set title "NS+AV" font ",13"
plot D."/sod/ref.dat" u 1:2 w l ls 1 t "Exact", D."/sod/ns_av.dat" u 1:2 w p ls 6 t "NS+AV"

# Row 2: Velocity
set title ""; set yrange [-0.05:1.05]
set ylabel "v" font ",13"
plot D."/sod/ref.dat" u 1:3 w l ls 1 notitle, D."/sod/sph.dat" u 1:3 w p ls 2 notitle
plot D."/sod/ref.dat" u 1:3 w l ls 1 notitle, D."/sod/voronoi.dat" u 1:3 w p ls 3 notitle
plot D."/sod/ref.dat" u 1:3 w l ls 1 notitle, D."/sod/price.dat" u 1:3 w p ls 4 notitle
plot D."/sod/ref.dat" u 1:3 w l ls 1 notitle, D."/sod/hllc_cd10.dat" u 1:3 w p ls 5 notitle
plot D."/sod/ref.dat" u 1:3 w l ls 1 notitle, D."/sod/ns_av.dat" u 1:3 w p ls 6 notitle

# Row 3: Pressure
set yrange [0.05:1.1]
set ylabel "P" font ",13"
plot D."/sod/ref.dat" u 1:4 w l ls 1 notitle, D."/sod/sph.dat" u 1:4 w p ls 2 notitle
plot D."/sod/ref.dat" u 1:4 w l ls 1 notitle, D."/sod/voronoi.dat" u 1:4 w p ls 3 notitle
plot D."/sod/ref.dat" u 1:4 w l ls 1 notitle, D."/sod/price.dat" u 1:4 w p ls 4 notitle
plot D."/sod/ref.dat" u 1:4 w l ls 1 notitle, D."/sod/hllc_cd10.dat" u 1:4 w p ls 5 notitle
plot D."/sod/ref.dat" u 1:4 w l ls 1 notitle, D."/sod/ns_av.dat" u 1:4 w p ls 6 notitle

# Row 4: Internal Energy
set yrange [1.5:3.2]
set xlabel "x" font ",12"
set ylabel "e" font ",13"
plot D."/sod/ref.dat" u 1:5 w l ls 1 notitle, D."/sod/sph.dat" u 1:5 w p ls 2 notitle
plot D."/sod/ref.dat" u 1:5 w l ls 1 notitle, D."/sod/voronoi.dat" u 1:5 w p ls 3 notitle
plot D."/sod/ref.dat" u 1:5 w l ls 1 notitle, D."/sod/price.dat" u 1:5 w p ls 4 notitle
plot D."/sod/ref.dat" u 1:5 w l ls 1 notitle, D."/sod/hllc_cd10.dat" u 1:5 w p ls 5 notitle
plot D."/sod/ref.dat" u 1:5 w l ls 1 notitle, D."/sod/ns_av.dat" u 1:5 w p ls 6 notitle
unset multiplot
set xlabel ""

# ============================================================
#  Interacting Blast Waves
# ============================================================
set output "1d_blast.pdf"
set multiplot layout 4,5 margins 0.07,0.97,0.06,0.95 spacing 0.045,0.04

# Row 1: Density
set xrange [0:1]; set yrange [0:13]
set ylabel "{/Symbol r}" font ",13"
set title "SPH" font ",13"
plot D."/blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/blast/sph.dat" u 1:2 w p ls 2 t "SPH"
set title "Voronoi" font ",13"
plot D."/blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/blast/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot D."/blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/blast/price.dat" u 1:2 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot D."/blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/blast/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"
set title "NS+AV" font ",13"
plot D."/blast/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/blast/ns_av.dat" u 1:2 w p ls 6 t "NS+AV"

# Row 2: Velocity
set title ""; set yrange [-1:16]
set ylabel "v" font ",13"
plot D."/blast/ref.dat" u 1:3 w l ls 1 notitle, D."/blast/sph.dat" u 1:3 w p ls 2 notitle
plot D."/blast/ref.dat" u 1:3 w l ls 1 notitle, D."/blast/voronoi.dat" u 1:3 w p ls 3 notitle
plot D."/blast/ref.dat" u 1:3 w l ls 1 notitle, D."/blast/price.dat" u 1:3 w p ls 4 notitle
plot D."/blast/ref.dat" u 1:3 w l ls 1 notitle, D."/blast/hllc_cd10.dat" u 1:3 w p ls 5 notitle
plot D."/blast/ref.dat" u 1:3 w l ls 1 notitle, D."/blast/ns_av.dat" u 1:3 w p ls 6 notitle

# Row 3: Pressure
set yrange [0:460]
set ylabel "P" font ",13"
plot D."/blast/ref.dat" u 1:4 w l ls 1 notitle, D."/blast/sph.dat" u 1:4 w p ls 2 notitle
plot D."/blast/ref.dat" u 1:4 w l ls 1 notitle, D."/blast/voronoi.dat" u 1:4 w p ls 3 notitle
plot D."/blast/ref.dat" u 1:4 w l ls 1 notitle, D."/blast/price.dat" u 1:4 w p ls 4 notitle
plot D."/blast/ref.dat" u 1:4 w l ls 1 notitle, D."/blast/hllc_cd10.dat" u 1:4 w p ls 5 notitle
plot D."/blast/ref.dat" u 1:4 w l ls 1 notitle, D."/blast/ns_av.dat" u 1:4 w p ls 6 notitle

# Row 4: Internal Energy
set yrange [0:1500]
set xlabel "x" font ",12"
set ylabel "e" font ",13"
plot D."/blast/ref.dat" u 1:5 w l ls 1 notitle, D."/blast/sph.dat" u 1:5 w p ls 2 notitle
plot D."/blast/ref.dat" u 1:5 w l ls 1 notitle, D."/blast/voronoi.dat" u 1:5 w p ls 3 notitle
plot D."/blast/ref.dat" u 1:5 w l ls 1 notitle, D."/blast/price.dat" u 1:5 w p ls 4 notitle
plot D."/blast/ref.dat" u 1:5 w l ls 1 notitle, D."/blast/hllc_cd10.dat" u 1:5 w p ls 5 notitle
plot D."/blast/ref.dat" u 1:5 w l ls 1 notitle, D."/blast/ns_av.dat" u 1:5 w p ls 6 notitle
unset multiplot
set xlabel ""

# ============================================================
#  Noh Strong Shock
# ============================================================
set output "1d_noh.pdf"
set multiplot layout 4,5 margins 0.07,0.97,0.06,0.95 spacing 0.045,0.04

# Row 1: Density
set xrange [0:0.45]; set yrange [0:5.5]
set ylabel "{/Symbol r}" font ",13"
set title "SPH" font ",13"
plot D."/noh/ref.dat" u 1:2 w l ls 1 t "Exact", D."/noh/sph.dat" u 1:2 w p ls 2 t "SPH"
set title "Voronoi" font ",13"
plot D."/noh/ref.dat" u 1:2 w l ls 1 t "Exact", D."/noh/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot D."/noh/ref.dat" u 1:2 w l ls 1 t "Exact", D."/noh/price.dat" u 1:2 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot D."/noh/ref.dat" u 1:2 w l ls 1 t "Exact", D."/noh/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"
set title "NS+AV" font ",13"
plot D."/noh/ref.dat" u 1:2 w l ls 1 t "Exact", D."/noh/ns_av.dat" u 1:2 w p ls 6 t "NS+AV"

# Row 2: Velocity
set title ""; set yrange [-1.1:0.15]
set ylabel "v" font ",13"
plot D."/noh/ref.dat" u 1:3 w l ls 1 notitle, D."/noh/sph.dat" u 1:3 w p ls 2 notitle
plot D."/noh/ref.dat" u 1:3 w l ls 1 notitle, D."/noh/voronoi.dat" u 1:3 w p ls 3 notitle
plot D."/noh/ref.dat" u 1:3 w l ls 1 notitle, D."/noh/price.dat" u 1:3 w p ls 4 notitle
plot D."/noh/ref.dat" u 1:3 w l ls 1 notitle, D."/noh/hllc_cd10.dat" u 1:3 w p ls 5 notitle
plot D."/noh/ref.dat" u 1:3 w l ls 1 notitle, D."/noh/ns_av.dat" u 1:3 w p ls 6 notitle

# Row 3: Pressure
set yrange [-0.05:2.0]
set ylabel "P" font ",13"
plot D."/noh/ref.dat" u 1:4 w l ls 1 notitle, D."/noh/sph.dat" u 1:4 w p ls 2 notitle
plot D."/noh/ref.dat" u 1:4 w l ls 1 notitle, D."/noh/voronoi.dat" u 1:4 w p ls 3 notitle
plot D."/noh/ref.dat" u 1:4 w l ls 1 notitle, D."/noh/price.dat" u 1:4 w p ls 4 notitle
plot D."/noh/ref.dat" u 1:4 w l ls 1 notitle, D."/noh/hllc_cd10.dat" u 1:4 w p ls 5 notitle
plot D."/noh/ref.dat" u 1:4 w l ls 1 notitle, D."/noh/ns_av.dat" u 1:4 w p ls 6 notitle

# Row 4: Internal Energy
set yrange [-0.02:1.1]
set xlabel "x" font ",12"
set ylabel "e" font ",13"
plot D."/noh/ref.dat" u 1:5 w l ls 1 notitle, D."/noh/sph.dat" u 1:5 w p ls 2 notitle
plot D."/noh/ref.dat" u 1:5 w l ls 1 notitle, D."/noh/voronoi.dat" u 1:5 w p ls 3 notitle
plot D."/noh/ref.dat" u 1:5 w l ls 1 notitle, D."/noh/price.dat" u 1:5 w p ls 4 notitle
plot D."/noh/ref.dat" u 1:5 w l ls 1 notitle, D."/noh/hllc_cd10.dat" u 1:5 w p ls 5 notitle
plot D."/noh/ref.dat" u 1:5 w l ls 1 notitle, D."/noh/ns_av.dat" u 1:5 w p ls 6 notitle
unset multiplot
set xlabel ""

# ============================================================
#  Shu-Osher Shock-Entropy Interaction
# ============================================================
set output "1d_shuosher.pdf"
set multiplot layout 4,5 margins 0.07,0.97,0.06,0.95 spacing 0.045,0.04

# Row 1: Density
set xrange [-5:5]; set yrange [0:4.8]
set ylabel "{/Symbol r}" font ",13"
set title "SPH" font ",13"
plot D."/shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/shuosher/sph.dat" u 1:2 w p ls 2 t "SPH"
set title "Voronoi" font ",13"
plot D."/shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/shuosher/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot D."/shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/shuosher/price.dat" u 1:2 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot D."/shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/shuosher/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"
set title "NS+AV" font ",13"
plot D."/shuosher/ref.dat" u 1:2 w l ls 1 t "Ref(4000)", D."/shuosher/ns_av.dat" u 1:2 w p ls 6 t "NS+AV"

# Row 2: Velocity
set title ""; set yrange [-4.5:3.5]
set ylabel "v" font ",13"
plot D."/shuosher/ref.dat" u 1:3 w l ls 1 notitle, D."/shuosher/sph.dat" u 1:3 w p ls 2 notitle
plot D."/shuosher/ref.dat" u 1:3 w l ls 1 notitle, D."/shuosher/voronoi.dat" u 1:3 w p ls 3 notitle
plot D."/shuosher/ref.dat" u 1:3 w l ls 1 notitle, D."/shuosher/price.dat" u 1:3 w p ls 4 notitle
plot D."/shuosher/ref.dat" u 1:3 w l ls 1 notitle, D."/shuosher/hllc_cd10.dat" u 1:3 w p ls 5 notitle
plot D."/shuosher/ref.dat" u 1:3 w l ls 1 notitle, D."/shuosher/ns_av.dat" u 1:3 w p ls 6 notitle

# Row 3: Pressure
set yrange [0:12]
set ylabel "P" font ",13"
plot D."/shuosher/ref.dat" u 1:4 w l ls 1 notitle, D."/shuosher/sph.dat" u 1:4 w p ls 2 notitle
plot D."/shuosher/ref.dat" u 1:4 w l ls 1 notitle, D."/shuosher/voronoi.dat" u 1:4 w p ls 3 notitle
plot D."/shuosher/ref.dat" u 1:4 w l ls 1 notitle, D."/shuosher/price.dat" u 1:4 w p ls 4 notitle
plot D."/shuosher/ref.dat" u 1:4 w l ls 1 notitle, D."/shuosher/hllc_cd10.dat" u 1:4 w p ls 5 notitle
plot D."/shuosher/ref.dat" u 1:4 w l ls 1 notitle, D."/shuosher/ns_av.dat" u 1:4 w p ls 6 notitle

# Row 4: Internal Energy
set yrange [0:13]
set xlabel "x" font ",12"
set ylabel "e" font ",13"
plot D."/shuosher/ref.dat" u 1:5 w l ls 1 notitle, D."/shuosher/sph.dat" u 1:5 w p ls 2 notitle
plot D."/shuosher/ref.dat" u 1:5 w l ls 1 notitle, D."/shuosher/voronoi.dat" u 1:5 w p ls 3 notitle
plot D."/shuosher/ref.dat" u 1:5 w l ls 1 notitle, D."/shuosher/price.dat" u 1:5 w p ls 4 notitle
plot D."/shuosher/ref.dat" u 1:5 w l ls 1 notitle, D."/shuosher/hllc_cd10.dat" u 1:5 w p ls 5 notitle
plot D."/shuosher/ref.dat" u 1:5 w l ls 1 notitle, D."/shuosher/ns_av.dat" u 1:5 w p ls 6 notitle
unset multiplot

# ============================================================
#  Lax Shock Tube
# ============================================================
set output "1d_lax.pdf"
set multiplot layout 4,5 margins 0.07,0.97,0.06,0.95 spacing 0.045,0.04

# Row 1: Density
set xrange [0:1]; set yrange [0.2:1.5]
set ylabel "{/Symbol r}" font ",13"
set title "SPH" font ",13"
plot D."/lax/ref.dat" u 1:2 w l ls 1 t "Exact", D."/lax/sph.dat" u 1:2 w p ls 2 t "SPH"
set title "Voronoi" font ",13"
plot D."/lax/ref.dat" u 1:2 w l ls 1 t "Exact", D."/lax/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot D."/lax/ref.dat" u 1:2 w l ls 1 t "Exact", D."/lax/price.dat" u 1:2 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot D."/lax/ref.dat" u 1:2 w l ls 1 t "Exact", D."/lax/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"
set title "NS+AV" font ",13"
plot D."/lax/ref.dat" u 1:2 w l ls 1 t "Exact", D."/lax/ns_av.dat" u 1:2 w p ls 6 t "NS+AV"

# Row 2: Velocity
set title ""; set yrange [-0.15:1.7]
set ylabel "v" font ",13"
plot D."/lax/ref.dat" u 1:3 w l ls 1 notitle, D."/lax/sph.dat" u 1:3 w p ls 2 notitle
plot D."/lax/ref.dat" u 1:3 w l ls 1 notitle, D."/lax/voronoi.dat" u 1:3 w p ls 3 notitle
plot D."/lax/ref.dat" u 1:3 w l ls 1 notitle, D."/lax/price.dat" u 1:3 w p ls 4 notitle
plot D."/lax/ref.dat" u 1:3 w l ls 1 notitle, D."/lax/hllc_cd10.dat" u 1:3 w p ls 5 notitle
plot D."/lax/ref.dat" u 1:3 w l ls 1 notitle, D."/lax/ns_av.dat" u 1:3 w p ls 6 notitle

# Row 3: Pressure
set yrange [0.3:4.0]
set ylabel "P" font ",13"
plot D."/lax/ref.dat" u 1:4 w l ls 1 notitle, D."/lax/sph.dat" u 1:4 w p ls 2 notitle
plot D."/lax/ref.dat" u 1:4 w l ls 1 notitle, D."/lax/voronoi.dat" u 1:4 w p ls 3 notitle
plot D."/lax/ref.dat" u 1:4 w l ls 1 notitle, D."/lax/price.dat" u 1:4 w p ls 4 notitle
plot D."/lax/ref.dat" u 1:4 w l ls 1 notitle, D."/lax/hllc_cd10.dat" u 1:4 w p ls 5 notitle
plot D."/lax/ref.dat" u 1:4 w l ls 1 notitle, D."/lax/ns_av.dat" u 1:4 w p ls 6 notitle

# Row 4: Internal Energy
set yrange [2:22]
set xlabel "x" font ",12"
set ylabel "e" font ",13"
plot D."/lax/ref.dat" u 1:5 w l ls 1 notitle, D."/lax/sph.dat" u 1:5 w p ls 2 notitle
plot D."/lax/ref.dat" u 1:5 w l ls 1 notitle, D."/lax/voronoi.dat" u 1:5 w p ls 3 notitle
plot D."/lax/ref.dat" u 1:5 w l ls 1 notitle, D."/lax/price.dat" u 1:5 w p ls 4 notitle
plot D."/lax/ref.dat" u 1:5 w l ls 1 notitle, D."/lax/hllc_cd10.dat" u 1:5 w p ls 5 notitle
plot D."/lax/ref.dat" u 1:5 w l ls 1 notitle, D."/lax/ns_av.dat" u 1:5 w p ls 6 notitle
unset multiplot
set xlabel ""

# ============================================================
#  Double Rarefaction (Toro Test #2)
# ============================================================
set output "1d_dblfan.pdf"
set multiplot layout 4,5 margins 0.07,0.97,0.06,0.95 spacing 0.045,0.04

# Row 1: Density
set xrange [0:1]; set yrange [0:1.15]
set ylabel "{/Symbol r}" font ",13"
set title "SPH" font ",13"
plot D."/dblfan/ref.dat" u 1:2 w l ls 1 t "Exact", D."/dblfan/sph.dat" u 1:2 w p ls 2 t "SPH"
set title "Voronoi" font ",13"
plot D."/dblfan/ref.dat" u 1:2 w l ls 1 t "Exact", D."/dblfan/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot D."/dblfan/ref.dat" u 1:2 w l ls 1 t "Exact", D."/dblfan/price.dat" u 1:2 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot D."/dblfan/ref.dat" u 1:2 w l ls 1 t "Exact", D."/dblfan/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"
set title "NS+AV" font ",13"
plot D."/dblfan/ref.dat" u 1:2 w l ls 1 t "Exact", D."/dblfan/ns_av.dat" u 1:2 w p ls 6 t "NS+AV"

# Row 2: Velocity
set title ""; set yrange [-2.3:2.3]
set ylabel "v" font ",13"
plot D."/dblfan/ref.dat" u 1:3 w l ls 1 notitle, D."/dblfan/sph.dat" u 1:3 w p ls 2 notitle
plot D."/dblfan/ref.dat" u 1:3 w l ls 1 notitle, D."/dblfan/voronoi.dat" u 1:3 w p ls 3 notitle
plot D."/dblfan/ref.dat" u 1:3 w l ls 1 notitle, D."/dblfan/price.dat" u 1:3 w p ls 4 notitle
plot D."/dblfan/ref.dat" u 1:3 w l ls 1 notitle, D."/dblfan/hllc_cd10.dat" u 1:3 w p ls 5 notitle
plot D."/dblfan/ref.dat" u 1:3 w l ls 1 notitle, D."/dblfan/ns_av.dat" u 1:3 w p ls 6 notitle

# Row 3: Pressure
set yrange [0:0.45]
set ylabel "P" font ",13"
plot D."/dblfan/ref.dat" u 1:4 w l ls 1 notitle, D."/dblfan/sph.dat" u 1:4 w p ls 2 notitle
plot D."/dblfan/ref.dat" u 1:4 w l ls 1 notitle, D."/dblfan/voronoi.dat" u 1:4 w p ls 3 notitle
plot D."/dblfan/ref.dat" u 1:4 w l ls 1 notitle, D."/dblfan/price.dat" u 1:4 w p ls 4 notitle
plot D."/dblfan/ref.dat" u 1:4 w l ls 1 notitle, D."/dblfan/hllc_cd10.dat" u 1:4 w p ls 5 notitle
plot D."/dblfan/ref.dat" u 1:4 w l ls 1 notitle, D."/dblfan/ns_av.dat" u 1:4 w p ls 6 notitle

# Row 4: Internal Energy
set yrange [0:1.15]
set xlabel "x" font ",12"
set ylabel "e" font ",13"
plot D."/dblfan/ref.dat" u 1:5 w l ls 1 notitle, D."/dblfan/sph.dat" u 1:5 w p ls 2 notitle
plot D."/dblfan/ref.dat" u 1:5 w l ls 1 notitle, D."/dblfan/voronoi.dat" u 1:5 w p ls 3 notitle
plot D."/dblfan/ref.dat" u 1:5 w l ls 1 notitle, D."/dblfan/price.dat" u 1:5 w p ls 4 notitle
plot D."/dblfan/ref.dat" u 1:5 w l ls 1 notitle, D."/dblfan/hllc_cd10.dat" u 1:5 w p ls 5 notitle
plot D."/dblfan/ref.dat" u 1:5 w l ls 1 notitle, D."/dblfan/ns_av.dat" u 1:5 w p ls 6 notitle
unset multiplot
set xlabel ""

# ============================================================
#  Two-Shock Collision (Toro Test #4)
# ============================================================
set output "1d_collision.pdf"
set multiplot layout 4,5 margins 0.07,0.97,0.06,0.95 spacing 0.045,0.04

# Row 1: Density
set xrange [0:1]; set yrange [3:34]
set ylabel "{/Symbol r}" font ",13"
set title "SPH" font ",13"
plot D."/collision/ref.dat" u 1:2 w l ls 1 t "Exact", D."/collision/sph.dat" u 1:2 w p ls 2 t "SPH"
set title "Voronoi" font ",13"
plot D."/collision/ref.dat" u 1:2 w l ls 1 t "Exact", D."/collision/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot D."/collision/ref.dat" u 1:2 w l ls 1 t "Exact", D."/collision/price.dat" u 1:2 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot D."/collision/ref.dat" u 1:2 w l ls 1 t "Exact", D."/collision/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"
set title "NS+AV" font ",13"
plot D."/collision/ref.dat" u 1:2 w l ls 1 t "Exact", D."/collision/ns_av.dat" u 1:2 w p ls 6 t "NS+AV"

# Row 2: Velocity
set title ""; set yrange [-8:22]
set ylabel "v" font ",13"
plot D."/collision/ref.dat" u 1:3 w l ls 1 notitle, D."/collision/sph.dat" u 1:3 w p ls 2 notitle
plot D."/collision/ref.dat" u 1:3 w l ls 1 notitle, D."/collision/voronoi.dat" u 1:3 w p ls 3 notitle
plot D."/collision/ref.dat" u 1:3 w l ls 1 notitle, D."/collision/price.dat" u 1:3 w p ls 4 notitle
plot D."/collision/ref.dat" u 1:3 w l ls 1 notitle, D."/collision/hllc_cd10.dat" u 1:3 w p ls 5 notitle
plot D."/collision/ref.dat" u 1:3 w l ls 1 notitle, D."/collision/ns_av.dat" u 1:3 w p ls 6 notitle

# Row 3: Pressure
set yrange [0:1850]
set ylabel "P" font ",13"
plot D."/collision/ref.dat" u 1:4 w l ls 1 notitle, D."/collision/sph.dat" u 1:4 w p ls 2 notitle
plot D."/collision/ref.dat" u 1:4 w l ls 1 notitle, D."/collision/voronoi.dat" u 1:4 w p ls 3 notitle
plot D."/collision/ref.dat" u 1:4 w l ls 1 notitle, D."/collision/price.dat" u 1:4 w p ls 4 notitle
plot D."/collision/ref.dat" u 1:4 w l ls 1 notitle, D."/collision/hllc_cd10.dat" u 1:4 w p ls 5 notitle
plot D."/collision/ref.dat" u 1:4 w l ls 1 notitle, D."/collision/ns_av.dat" u 1:4 w p ls 6 notitle

# Row 4: Internal Energy
set yrange [0:400]
set xlabel "x" font ",12"
set ylabel "e" font ",13"
plot D."/collision/ref.dat" u 1:5 w l ls 1 notitle, D."/collision/sph.dat" u 1:5 w p ls 2 notitle
plot D."/collision/ref.dat" u 1:5 w l ls 1 notitle, D."/collision/voronoi.dat" u 1:5 w p ls 3 notitle
plot D."/collision/ref.dat" u 1:5 w l ls 1 notitle, D."/collision/price.dat" u 1:5 w p ls 4 notitle
plot D."/collision/ref.dat" u 1:5 w l ls 1 notitle, D."/collision/hllc_cd10.dat" u 1:5 w p ls 5 notitle
plot D."/collision/ref.dat" u 1:5 w l ls 1 notitle, D."/collision/ns_av.dat" u 1:5 w p ls 6 notitle
unset multiplot
set xlabel ""

# ============================================================
#  Stationary Contact Discontinuity
# ============================================================
set output "1d_contact.pdf"
set multiplot layout 4,5 margins 0.07,0.97,0.06,0.95 spacing 0.045,0.04

# Row 1: Density
set xrange [0:1]; set yrange [0.05:1.6]
set ylabel "{/Symbol r}" font ",13"
set title "SPH" font ",13"
plot D."/contact/ref.dat" u 1:2 w l ls 1 t "Exact", D."/contact/sph.dat" u 1:2 w p ls 2 t "SPH"
set title "Voronoi" font ",13"
plot D."/contact/ref.dat" u 1:2 w l ls 1 t "Exact", D."/contact/voronoi.dat" u 1:2 w p ls 3 t "Voronoi"
set title "Price AV" font ",13"
plot D."/contact/ref.dat" u 1:2 w l ls 1 t "Exact", D."/contact/price.dat" u 1:2 w p ls 4 t "Price AV"
set title "HLLC+CD10" font ",13"
plot D."/contact/ref.dat" u 1:2 w l ls 1 t "Exact", D."/contact/hllc_cd10.dat" u 1:2 w p ls 5 t "HLLC+CD10"
set title "NS+AV" font ",13"
plot D."/contact/ref.dat" u 1:2 w l ls 1 t "Exact", D."/contact/ns_av.dat" u 1:2 w p ls 6 t "NS+AV"

# Row 2: Velocity
set title ""; set yrange [-0.06:0.06]
set ylabel "v" font ",13"
plot D."/contact/ref.dat" u 1:3 w l ls 1 notitle, D."/contact/sph.dat" u 1:3 w p ls 2 notitle
plot D."/contact/ref.dat" u 1:3 w l ls 1 notitle, D."/contact/voronoi.dat" u 1:3 w p ls 3 notitle
plot D."/contact/ref.dat" u 1:3 w l ls 1 notitle, D."/contact/price.dat" u 1:3 w p ls 4 notitle
plot D."/contact/ref.dat" u 1:3 w l ls 1 notitle, D."/contact/hllc_cd10.dat" u 1:3 w p ls 5 notitle
plot D."/contact/ref.dat" u 1:3 w l ls 1 notitle, D."/contact/ns_av.dat" u 1:3 w p ls 6 notitle

# Row 3: Pressure
set yrange [0.85:1.06]
set ylabel "P" font ",13"
plot D."/contact/ref.dat" u 1:4 w l ls 1 notitle, D."/contact/sph.dat" u 1:4 w p ls 2 notitle
plot D."/contact/ref.dat" u 1:4 w l ls 1 notitle, D."/contact/voronoi.dat" u 1:4 w p ls 3 notitle
plot D."/contact/ref.dat" u 1:4 w l ls 1 notitle, D."/contact/price.dat" u 1:4 w p ls 4 notitle
plot D."/contact/ref.dat" u 1:4 w l ls 1 notitle, D."/contact/hllc_cd10.dat" u 1:4 w p ls 5 notitle
plot D."/contact/ref.dat" u 1:4 w l ls 1 notitle, D."/contact/ns_av.dat" u 1:4 w p ls 6 notitle

# Row 4: Internal Energy
set yrange [1.5:22]
set xlabel "x" font ",12"
set ylabel "e" font ",13"
plot D."/contact/ref.dat" u 1:5 w l ls 1 notitle, D."/contact/sph.dat" u 1:5 w p ls 2 notitle
plot D."/contact/ref.dat" u 1:5 w l ls 1 notitle, D."/contact/voronoi.dat" u 1:5 w p ls 3 notitle
plot D."/contact/ref.dat" u 1:5 w l ls 1 notitle, D."/contact/price.dat" u 1:5 w p ls 4 notitle
plot D."/contact/ref.dat" u 1:5 w l ls 1 notitle, D."/contact/hllc_cd10.dat" u 1:5 w p ls 5 notitle
plot D."/contact/ref.dat" u 1:5 w l ls 1 notitle, D."/contact/ns_av.dat" u 1:5 w p ls 6 notitle
unset multiplot
