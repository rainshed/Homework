set xlabel 'X'
set ylabel 'Wave function'
plot 'data/p3_diag.dat'
set term pngcairo
set output "picture/p3_diag.png"
replot
set output
