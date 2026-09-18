# plot_Cd_body0.gp  (run from case root)
# Input: postProcessing/forces/dragFx_body0.dat  columns: time Fx

FILE = "postProcessing/forces/dragFx_body0.dat"

rho  = 1000.0
Uinf = 1.0
Aref = 0.1
den  = 0.5*rho*Uinf*Uinf*Aref   # = 50

set term pngcairo size 1200,600
set output "postProcessing/forces/Cd_body0.png"

set title "Cd vs time (body0)"
set xlabel "time"
set ylabel "Cd = Fx / (0.5*rho*U^2*Aref)"
set grid
set key left top
set datafile commentschars "#"

plot FILE using 1:($2/den) with lines lw 2 title sprintf("den=%.3g", den)

set output
print "Wrote postProcessing/forces/Cd_body0.png"
