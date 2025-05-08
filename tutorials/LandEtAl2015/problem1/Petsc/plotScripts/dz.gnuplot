set term pdfcairo dashed enhanced
set datafile separator " "
set output "dz_vs_cellSize.pdf"

set size ratio -1
set grid

# Manual axis ranges
set xrange [0.5:4.5]
set yrange [2:5]

# xtics are labels, indexed by position (1 to 4)
set xtics ("0.5" 1, "0.25" 2, "0.125" 3, "0.0625" 4)

set xlabel "Cell size (in mm)"
set ylabel "Dz (in mm)"
set key outside left center

# Plot using index as x (1 to 4), dz (scaled to mm) as y
plot "dz_values.dat" using ($0+1):(1e3*$2) with linespoints pt 6 ps 1 lc "blue" title "Dz displacement"
