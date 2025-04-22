set term pdfcairo dashed enhanced
set datafile separator " "

set output "poly.midline.pdf"

set size ratio -1

set grid
set xrange [-8.75:0]
set yrange [-19:5]
set xtics 2
#set xtics add (5, 25, 50)
set ytics
#set logscale x
#set logscale y
#set format y "10^{%L}"
#set ytics 0.002
set xlabel "x (in mm)"
set ylabel "z (in mm)"
#set key left top;
#set key right bottom;
set key outside left center;

#set label "1^{st} order" at graph 0.5,0.86 center rotate by 10
#set label "2^{nd} order" at graph 0.5,0.37 center rotate by 25

# Calculate the deltaX as cbrt(totalVolume/numCells)
plot \
    "poly.hypre.1/midLineDeformed.txt" u (1e3*$1):(1e3*$3) every 9::0 w lp pt 6 ps 1 lc "blue" t "{/Symbol \D}x = 1.26 mm", \
    "poly.hypre.2/midLineDeformed.txt" u (1e3*$1):(1e3*$3) every 9::0 w lp pt 6 ps 0.8 lc "cyan" t "{/Symbol \D}x = 0.63 mm", \
    "poly.hypre.3/midLineDeformed.txt" u (1e3*$1):(1e3*$3) every 9::0 w lp pt 6 ps 0.6 lc "green" t "{/Symbol \D}x = 0.32 mm", \
    "poly.hypre.4/midLineDeformed.txt" u (1e3*$1):(1e3*$3) every 9::0 w lp pt 6 ps 0.4 lc "red" t "{/Symbol \D}x = 0.16 mm"

# Apex plot
set output "midline_apex.pdf"
set xrange [-2:0]
set yrange [-15:-13]
set xtics 2
set size ratio 1
plot \
    "poly.hypre.1/midLineDeformed.txt" u (1e3*$1):(1e3*$3) w lp pt 6 ps 1 lc "blue" t "{/Symbol \D}x = 1.26 mm", \
    "poly.hypre.2/midLineDeformed.txt" u (1e3*$1):(1e3*$3) w lp pt 6 ps 0.8 lc "cyan" t "{/Symbol \D}x = 0.63 mm", \
    "poly.hypre.3/midLineDeformed.txt" u (1e3*$1):(1e3*$3) w lp pt 6 ps 0.6 lc "green" t "{/Symbol \D}x = 0.32 mm", \
    "poly.hypre.4/midLineDeformed.txt" u (1e3*$1):(1e3*$3) w lp pt 6 ps 0.4 lc "red" t "{/Symbol \D}x = 0.16 mm"

# Inflection plot
set output "midline_inflection.pdf"
set xrange [-8.75:-8.25]
set yrange [-2:2]
set xtics 0.25
set size ratio 1
plot \
    "poly.hypre.1/midLineDeformed.txt" u (1e3*$1):(1e3*$3) w lp pt 6 ps 1 lc "blue" t "{/Symbol \D}x = 1.26 mm", \
    "poly.hypre.2/midLineDeformed.txt" u (1e3*$1):(1e3*$3) w lp pt 6 ps 0.8 lc "cyan" t "{/Symbol \D}x = 0.63 mm", \
    "poly.hypre.3/midLineDeformed.txt" u (1e3*$1):(1e3*$3) w lp pt 6 ps 0.6 lc "green" t "{/Symbol \D}x = 0.32 mm", \
    "poly.hypre.4/midLineDeformed.txt" u (1e3*$1):(1e3*$3) w lp pt 6 ps 0.4 lc "red" t "{/Symbol \D}x = 0.16 mm"
