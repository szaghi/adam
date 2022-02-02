set terminal png enhanced notransparent rounded large size 1600,1200 nocrop linewidth 2 background "#ffffff" font "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf,20"
set output 'BNAME.png'
set autoscale
set yrange [1:15]
set title 'Slice along sphere center-axis'
set xlabel 'X'
set ylabel 'rho'
set border 3
set tics nomirror 
set grid
plot "BNAME.dat" using 3:6:1 ev 8 with labels offset 0,1 notitle, "" using 3:6 with linespoints lc 7 lw 2 pt 7 notitle
