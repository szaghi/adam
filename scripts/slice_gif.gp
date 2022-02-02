set terminal gif enhanced notransparent rounded large size 1600,1200 nocrop linewidth 2 background "#ffffff" font "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf,20" animate delay 10
set output 'slice.gif'
set autoscale
set yrange [1:12]
set title 'Slice along sphere center-axis'
set xlabel 'X'
set ylabel 'rho'
set border 3
set tics nomirror 
set grid
do for [fn in system("ls *slice_*dat")] {plot fn using 3:6:1 ev 8 with labels offset 0,1 notitle, "" using 3:6 with linespoints lc 7 lw 2 pt 7 notitle; pause 0.5 }
