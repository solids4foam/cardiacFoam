set datafile commentschars "#"

set terminal pdfcairo enhanced font "Times,18" size 6in,4in
set output "Cd_vs_time.pdf"

set title "Drag Coefficient vs Time"
set xlabel "Time"
set ylabel "Cd"

set grid
set border lw 1.5
set key top left

plot "postProcessing/forceCoeffs/0/coefficient.dat" using 1:2 with lines lw 2 title "Cd"