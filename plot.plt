# plot.plt
set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*"
set title "Fortran Example"
set grid
set xlabel "x"
set ylabel "y"
m="test_data.txt"
m_prime="test_data_prime.txt"
plot m using 1:2 with linespoints, m_prime using 1:2 with linespoints
