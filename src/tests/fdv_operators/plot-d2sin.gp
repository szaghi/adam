set terminal wxt size 1200,800

set multiplot layout 4,4 title "Derivative 2, numerical FDV vs with exact solution" font ",14"

set title "FD CC O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_02-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FD_CC_order_02-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FD CC O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_04-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FD_CC_order_04-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FD CC O 06"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_06-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FD_CC_order_06-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FD CC O 08"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_08-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FD_CC_order_08-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FD CC O 10"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_10-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FD_CC_order_10-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV CC O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_02-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FV_CC_order_02-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV CC O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_04-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FV_CC_order_04-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV CC O 06"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_06-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FV_CC_order_06-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV CC O 08"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_08-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FV_CC_order_08-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV LU O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_LU_order_02-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FV_LU_order_02-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV LU O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_LU_order_04-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FV_LU_order_04-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV RU O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_RU_order_02-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FV_RU_order_02-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

set title "FV RU O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_RU_order_04-d2sin.dat' skip 1 using 1:2 with lines lw 2 title 'd2sin', 'FV_RU_order_04-d2sin.dat' skip 1 using 1:3 with lines lw 2 title 'cos'

unset multiplot

# save png
set terminal png size 2400,1600
set output 'fdv-d2sin.png'
replot
