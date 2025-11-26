set terminal wxt size 1200,800

set multiplot layout 4,4 title "Derivative 1, numerical FDV vs with exact solution" font ",14"

set title "FD CC O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_02-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FD_CC_order_02-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FD CC O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_04-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FD_CC_order_04-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FD CC O 06"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_06-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FD_CC_order_06-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FD CC O 08"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_08-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FD_CC_order_08-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FD CC O 10"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FD_CC_order_10-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FD_CC_order_10-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV CC O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_02-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_CC_order_02-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV CC O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_04-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_CC_order_04-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV CC O 06"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_06-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_CC_order_06-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV CC O 08"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_08-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_CC_order_08-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV CC O 10"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_CC_order_10-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_CC_order_10-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV LU O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_LU_order_02-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_LU_order_02-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV LU O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_LU_order_04-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_LU_order_04-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV RU O 02"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_RU_order_02-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_RU_order_02-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

set title "FV RU O 04"
set xlabel "x"
set ylabel "y"
set grid
set key top right
plot 'FV_RU_order_04-step.dat' skip 1 using 1:2 with lines lw 2 title 'numerical', 'FV_RU_order_04-step.dat' skip 1 using 1:3 with lines lw 2 title 'exact'

unset multiplot

# save png
set terminal png size 2400,1600
set output 'fdv-step.png'
replot
