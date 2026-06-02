set terminal pdfcairo enhanced color font "Arial,18" size 10,6
set output "Cd_comparison.pdf"

set datafile commentschars "#"

set title "Nonlinear bending of beam" font "Arial-Bold,22"

set xlabel "Time [s]" font "Arial,20"
set ylabel "Cd" font "Arial,20"

set xrange [0:8]
set yrange [-80:100]

set grid
set border linewidth 1.5
set tics font "Arial,18"
set key inside top right box opaque font "Arial,16" spacing 1.3

plot \
"embedded/coefficient.dat" using 1:2 with lines lw 3 lc rgb "black" title "Embedded", \
"immersed/Cd_body0.dat" using 1:2 with lines lw 3 lc rgb "blue" title "Immersed", \
"immersedMeshMotion/Cd_body0.dat" using 1:2 with lines lw 3 lc rgb "red" title "Immersed + Mesh Motion"