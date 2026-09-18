set terminal pdfcairo size 12in,7in font ",14"
set output "Cd_comparison_400k_3p2M_reference.pdf"

set grid
set xlabel "t [s]"
set ylabel "C_d [-]"
set title "Oscillating sphere: drag coefficient comparison"

set key right top
set border lw 1.5
set datafile commentschars "#"

set xrange [0:10]
set yrange [-30:30]

plot \
  "Cd_body0_400k.dat" using 1:2 with lines lw 2 lc rgb "blue" title "Present: 400k cells", \
  "Cd_body0_800k.dat" using 1:2 with lines lw 2 lc rgb "black" title "Present: 3.2M cells", \
  "reference_Gilmanov_Sotiropoulos_Cd_sorted.dat" using 1:2 with points pt 7 ps 0.55 lc rgb "red" title "Gilmanov \& Sotiropoulos"

set output