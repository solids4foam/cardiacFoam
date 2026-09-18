reset

set datafile commentschars "#"

set terminal pdfcairo enhanced font "Times,18" size 6in,4in
set output "Cd_IB_coefficient_comparison.pdf"

set title "Oscillating Cylinder Drag Coefficient Comparison"
set xlabel "Time"
set ylabel "C_d"

set grid
set key top right
set border lw 1.5
set xrange [0:12]

# Same mesh, different IB coefficients: solid coloured lines
set style line 1 lc rgb "#4169E1" lw 2.5 dt 1   # royal blue
set style line 2 lc rgb "#2ca02c" lw 2.5 dt 1   # green
set style line 3 lc rgb "#000000" lw 2.5 dt 1   # black

# Reference: red dots
set style line 4 lc rgb "red" pt 7 ps 0.40

plot \
"Cd_body0-Mesh56320-IBPt05.dat" using 1:2 with lines  ls 1 title "IB coefficient = 0.05", \
"Cd_body0-Mesh56320-IBPt2.dat"  using 1:2 with lines  ls 2 title "IB coefficient = 0.2", \
"Cd_body0-Mesh56320-IBPt8.dat"  using 1:2 with lines  ls 3 title "IB coefficient = 0.8", \
"reference_Cd_sorted.dat"       using 1:2 with points ls 4 title "Wan \\& Turek"