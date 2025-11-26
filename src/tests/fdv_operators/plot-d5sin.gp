set terminal wxt size 1200,800

set multiplot layout 2,2 title "Derivative 5, numerical FDV vs with exact solution" font ",14"

set title "FD CC O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_02-d5sin.dat' skip 1 using 1:2 with lines lw 2 title 'd5sin', 'FD_CC_order_02-d5sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FD CC O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_04-d5sin.dat' skip 1 using 1:2 with lines lw 2 title 'd5sin', 'FD_CC_order_04-d5sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FD CC O 06"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_06-d5sin.dat' skip 1 using 1:2 with lines lw 2 title 'd5sin', 'FD_CC_order_06-d5sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV CC O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_02-d5sin.dat' skip 1 using 1:2 with lines lw 2 title 'd5sin', 'FV_CC_order_02-d5sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

unset multiplot

# save png
set terminal png size 2400,1600
set output 'fdv-d5sin.png'
replot
