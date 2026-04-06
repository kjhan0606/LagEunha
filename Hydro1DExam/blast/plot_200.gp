set terminal pngcairo size 2250,1200 enhanced font 'Arial,11'
set output 'blast/result_200.png'
set multiplot layout 4,5 title 'Interacting Blast Waves (N=200)' font ',14'
set key top right font ',9'

set title 'HLLC+CD10' font ',12'
set ylabel 'Density'
set xlabel 'x'
set xrange [0:1]
set yrange [0:13.4355]
plot 'blast/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/hllc_cd10.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#a65628' t 'HLLC+CD10'

set title 'Price AV' font ',12'
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:13.4355]
plot 'blast/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/price.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#e41a1c' t 'Price AV'

set title 'Vor M(3,4)+dpZ' font ',12'
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:13.4355]
plot 'blast/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/voronoi.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#4daf4a' t 'Vor M(3,4)+dpZ'

set title 'SPH' font ',12'
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:13.4355]
plot 'blast/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/sph.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#984ea3' t 'SPH'

set title 'NS+AV' font ',12'
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:13.4355]
plot 'blast/ref.dat' u 1:2 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/ns_av.dat' u 1:2 w p pt 2 ps 1.0 lw 2 lc rgb '#ff7f00' t 'NS+AV'

unset title
set ylabel 'Velocity'
set xlabel 'x'
set xrange [0:1]
set yrange [-1.2439:16.2283]
plot 'blast/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/hllc_cd10.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#a65628' t 'HLLC+CD10'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [-1.2439:16.2283]
plot 'blast/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/price.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#e41a1c' t 'Price AV'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [-1.2439:16.2283]
plot 'blast/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/voronoi.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#4daf4a' t 'Vor M(3,4)+dpZ'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [-1.2439:16.2283]
plot 'blast/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/sph.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#984ea3' t 'SPH'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [-1.2439:16.2283]
plot 'blast/ref.dat' u 1:3 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/ns_av.dat' u 1:3 w p pt 2 ps 1.0 lw 2 lc rgb '#ff7f00' t 'NS+AV'

unset title
set ylabel 'Pressure'
set xlabel 'x'
set xrange [0:1]
set yrange [0:488.476]
plot 'blast/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/hllc_cd10.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#a65628' t 'HLLC+CD10'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:488.476]
plot 'blast/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/price.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#e41a1c' t 'Price AV'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:488.476]
plot 'blast/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/voronoi.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#4daf4a' t 'Vor M(3,4)+dpZ'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:488.476]
plot 'blast/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/sph.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#984ea3' t 'SPH'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:488.476]
plot 'blast/ref.dat' u 1:4 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/ns_av.dat' u 1:4 w p pt 2 ps 1.0 lw 2 lc rgb '#ff7f00' t 'NS+AV'

unset title
set ylabel 'Internal Energy'
set xlabel 'x'
set xrange [0:1]
set yrange [0:707.205]
plot 'blast/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/hllc_cd10.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#a65628' t 'HLLC+CD10'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:707.205]
plot 'blast/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/price.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#e41a1c' t 'Price AV'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:707.205]
plot 'blast/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/voronoi.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#4daf4a' t 'Vor M(3,4)+dpZ'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:707.205]
plot 'blast/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/sph.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#984ea3' t 'SPH'

unset title
unset ylabel
set xlabel 'x'
set xrange [0:1]
set yrange [0:707.205]
plot 'blast/ref.dat' u 1:5 w l lw 2 lc rgb '#000000' t 'Ref', \
     'blast/ns_av.dat' u 1:5 w p pt 2 ps 1.0 lw 2 lc rgb '#ff7f00' t 'NS+AV'

unset multiplot
