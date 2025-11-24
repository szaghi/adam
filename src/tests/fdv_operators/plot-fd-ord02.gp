set terminal wxt size 1200,800

set multiplot layout 2,3 title "Numerical FD (order 2) vs with exact functions" font ",14"

set title "sin1 vs cos"
set xlabel "x"
set ylabel "y"
set grid
set key top right

plot 'fd_d1sin-ord02.dat' using 1:2 with lines lw 2 title 'sin1', \
     'fd_d1sin-ord02.dat' using 1:3 with lines lw 2 title 'cos'

set title "sin2 vs -sin"
set xlabel "x"
set ylabel "y"
set grid
set key top right

plot 'fd_d2sin-ord02.dat' using 1:2 with lines lw 2 title 'sin2', \
     'fd_d2sin-ord02.dat' using 1:3 with lines lw 2 title '-sin'

set title "sin3 vs -cos"
set xlabel "x"
set ylabel "y"
set grid
set key top right

plot 'fd_d3sin-ord02.dat' using 1:2 with lines lw 2 title 'sin3', \
     'fd_d3sin-ord02.dat' using 1:3 with lines lw 2 title '-cos'

set title "sin4 vs sin"
set xlabel "x"
set ylabel "y"
set grid
set key top right

plot 'fd_d4sin-ord02.dat' using 1:2 with lines lw 2 title 'sin4', \
     'fd_d4sin-ord02.dat' using 1:3 with lines lw 2 title 'sin'

set title "sin5 vs cos"
set xlabel "x"
set ylabel "y"
set grid
set key top right

plot 'fd_d5sin-ord02.dat' using 1:2 with lines lw 2 title 'sin5', \
     'fd_d5sin-ord02.dat' using 1:3 with lines lw 2 title 'cos'

set title "sin6 vs -sin"
set xlabel "x"
set ylabel "y"
set grid
set key top right

plot 'fd_d6sin-ord02.dat' using 1:2 with lines lw 2 title 'sin6', \
     'fd_d6sin-ord02.dat' using 1:3 with lines lw 2 title '-sin'

unset multiplot

# save png
set terminal png size 1200,800
set output 'fd-ord02.png'
replot
