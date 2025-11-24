set terminal wxt size 1200,800

set multiplot layout 1,2 title "Numerical FD (order 10) vs with exact functions" font ",14"

set title "sin1 vs cos"
set xlabel "x"
set ylabel "y"
set grid
set key top right

plot 'fd_d1sin-ord10.dat' using 1:2 with lines lw 2 title 'sin1', \
     'fd_d1sin-ord10.dat' using 1:3 with lines lw 2 title 'cos'

set title "sin2 vs -sin"
set xlabel "x"
set ylabel "y"
set grid
set key top right

plot 'fd_d2sin-ord10.dat' using 1:2 with lines lw 2 title 'sin2', \
     'fd_d2sin-ord10.dat' using 1:3 with lines lw 2 title '-sin'

unset multiplot

# save png
set terminal png size 1200,800
set output 'fd-ord10.png'
replot
