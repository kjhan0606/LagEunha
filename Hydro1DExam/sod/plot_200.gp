set terminal pngcairo size 2250,1200 enhanced font 'Arial,11'
set output 'sod/result_200.png'
set multiplot layout 4,5 title 'Sod Shock Tube (N=200)' font ',14'
set key top right font ',9'

set title 'HLLC+CD10' font ',12'
set ylabel 'Density'
set xlabel 'x'
set xrange [0:1]
set yrange [0.044168:1.07975]
plot 'sod/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/hllc_cd10.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#a65628' t 'HLLC+CD10'

set title 'HLLC+Price' font ',12'
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0.044168:1.07975]
plot 'sod/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/hllc_price.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#377eb8' t 'HLLC+Price'

set title 'Price AV' font ',12'
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0.044168:1.07975]
plot 'sod/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/price.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#e41a1c' t 'Price AV'

set title 'Vor M(3,4)+dpZ' font ',12'
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0.044168:1.07975]
plot 'sod/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/voronoi.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#4daf4a' t 'Vor M(3,4)+dpZ'

set title 'SPH' font ',12'
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0.044168:1.07975]
plot 'sod/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/sph.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#984ea3' t 'SPH'

unset title
set ylabel 'Velocity'
set xlabel 'x'
set xrange [0:1]
set yrange [-0.0907105:1.08275]
plot 'sod/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/hllc_cd10.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#a65628' t 'HLLC+CD10'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [-0.0907105:1.08275]
plot 'sod/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/hllc_price.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#377eb8' t 'HLLC+Price'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [-0.0907105:1.08275]
plot 'sod/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/price.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#e41a1c' t 'Price AV'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [-0.0907105:1.08275]
plot 'sod/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/voronoi.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#4daf4a' t 'Vor M(3,4)+dpZ'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [-0.0907105:1.08275]
plot 'sod/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/sph.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#984ea3' t 'SPH'

unset title
set ylabel 'Pressure'
set xlabel 'x'
set xrange [0:1]
set yrange [0.0266245:1.08443]
plot 'sod/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/hllc_cd10.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#a65628' t 'HLLC+CD10'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0.0266245:1.08443]
plot 'sod/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/hllc_price.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#377eb8' t 'HLLC+Price'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0.0266245:1.08443]
plot 'sod/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/price.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#e41a1c' t 'Price AV'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0.0266245:1.08443]
plot 'sod/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/voronoi.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#4daf4a' t 'Vor M(3,4)+dpZ'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0.0266245:1.08443]
plot 'sod/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/sph.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#984ea3' t 'SPH'

unset title
set ylabel 'Internal Energy'
set xlabel 'x'
set xrange [0:1]
set yrange [1.58106:3.14553]
plot 'sod/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/hllc_cd10.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#a65628' t 'HLLC+CD10'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [1.58106:3.14553]
plot 'sod/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/hllc_price.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#377eb8' t 'HLLC+Price'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [1.58106:3.14553]
plot 'sod/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/price.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#e41a1c' t 'Price AV'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [1.58106:3.14553]
plot 'sod/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/voronoi.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#4daf4a' t 'Vor M(3,4)+dpZ'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [1.58106:3.14553]
plot 'sod/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'sod/sph.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#984ea3' t 'SPH'

unset multiplot
