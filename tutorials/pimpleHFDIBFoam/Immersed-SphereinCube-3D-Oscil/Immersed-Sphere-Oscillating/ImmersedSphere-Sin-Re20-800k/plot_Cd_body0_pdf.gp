FILE = "postProcessing/forces/dragFx_body0.dat"

rho  = 1.0
f    = 0.25
A    = 0.25
Uinf = 2.0*pi*f*A          # 0.392699...
D    = 0.1
Lz   = 0.1                 # <-- set to your actual mesh thickness
Aref = D*Lz
den  = 0.5*rho*Uinf*Uinf*Aref

set term pdfcairo size 16cm,9cm font "Helvetica,10"
set output "postProcessing/forces/Cd_body0.pdf"

set title "Cd vs time (body0)"
set xlabel "t [s]"
set ylabel "Cd = Fx / (0.5*rho*Uinf^2*Aref)"
set grid
set key left top
set datafile commentschars "#"

plot FILE using 1:($2/den) with lines lw 2 title sprintf("den=%.6g", den)

set output
print "Wrote postProcessing/forces/Cd_body0.pdf"