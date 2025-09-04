set term pdfcairo dashed enhanced
set datafile separator " "
set output "dz_vs_cellSize_Bulk_Anslysis.pdf"

set size ratio -1
set grid

set xrange [0.5:4.5]
set yrange [0:5]
set xtics ("0.5" 1, "0.25" 2, "0.125" 3, "0.0625" 4)

set xlabel "Cell size (in mm)"
set ylabel "Dz (in mm)"
set key outside left center

# Define a dashed line style for the reference
set style line 99 lt 2 lc rgb "black" lw 2 dashtype 2

# Plot all 2e*.dat files with their names as titles
plot \
  "2e3.dat" using ($0+1):(1e3*$2) with linespoints pt 6 ps 1 lc "red" title "2e3", \
  "2e4.dat" using ($0+1):(1e3*$2) with linespoints pt 6 ps 1 lc "blue" title "2e4", \
  "2e5.dat" using ($0+1):(1e3*$2) with linespoints pt 6 ps 1 lc "green" title "2e5", \
  "2e6.dat" using ($0+1):(1e3*$2) with linespoints pt 6 ps 1 lc "purple" title "2e6", \
  '+' using 1:(4.18) with lines linestyle 99 title "Reference (4.18 mm)"
