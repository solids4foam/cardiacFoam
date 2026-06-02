# Plots the drag coefficient based on the data calculate by the `forceCoeffs`
# function object.
# ---------------------------------------------------------------------------- #
set terminal pdfcairo enhanced color solid

set datafile separator whitespace
set datafile commentschars "#"

if (ARGC != 3) {
    print "Error: Wrong number of input arguments."
    print "usage:", ARG0, "<datafile> <verificationDataFile> <outfile>"
    print "       dataFile      Path to the post-processing data"
    print "       CdDataFile    Path to the verification data"
    print "       outFile       Path to the output file"
    exit 1
} else {
    dataFile = ARG1
    CdDataFile = ARG2
    set output ARG3
}

set grid
set tics nomirror
set border 3
set key top left
set key box width 2 height 1

lineColor = "blue"
pointColor = "red"
pointType = 6 # Hollow circle
pointSize = .2

# --- Plot Drag Coefficient -------------------------------------------------- #

set title "Drag Coefficient (C_d) vs. Time (t)"

set xlabel "t [s]"
set ylabel "C_d [1]"

set xtics 1
set ytics 10

set xrange [0:10]
set yrange [-30:30]

plot \
    dataFile \
        using 1:2 \
        title "Present" \
        with lines lc rgb lineColor, \
    CdDataFile \
        using 1:2 \
        title "Gilmanov\\&Sotiropoulos" \
        with points pt pointType ps pointSize lc rgb pointColor
