# plot.plt
set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*"
set title "Fortran Example"
set grid
set xlabel "x"
set ylabel "y"
m="stepsizecorr_cs.txt"
plot m using 1:2 with linespoints
