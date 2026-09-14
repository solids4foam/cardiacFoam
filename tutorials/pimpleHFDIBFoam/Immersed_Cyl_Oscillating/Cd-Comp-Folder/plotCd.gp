reset

set datafile commentschars "#"

set terminal pdfcairo enhanced font "Times,18" size 6in,4in
set output "Cd_mesh_comparison_t12.pdf"

set title "Oscillating Cylinder Drag Coefficient Comparison"
set xlabel "Time"
set ylabel "C_d"

set grid
set key top right
set border lw 1.5
set xrange [0:12]

# 14080: royal blue solid
set style line 1 lc rgb "#4169E1" lw 2.5 dt 1

# 56320: royal blue dash-dotted
set style line 2 lc rgb "#4169E1" lw 2.5 dashtype (8,4,2,4)

# 225280: black dash-dotted
set style line 3 lc rgb "#000000" lw 2.5 dashtype (8,4,2,4)

# 901120: black solid
set style line 4 lc rgb "#000000" lw 2.5 dt 1

# Reference: red dots
set style line 5 lc rgb "red" pt 7 ps 0.40

plot \
"Cd_body0-Mesh14080.dat"  using 1:2 with lines  ls 1 title "Mesh 14080", \
"Cd_body0-Mesh56320.dat"  using 1:2 with lines  ls 2 title "Mesh 56320", \
"Cd_body0-Mesh225280.dat" using 1:2 with lines  ls 3 title "Mesh 225280", \
"Cd_body0-Mesh901120.dat" using 1:2 with lines  ls 4 title "Mesh 901120", \
"reference_Cd_sorted.dat" using 1:2 with points ls 5 title "Wan \\& Turek"