set xlabel 'X'
set ylabel 'Wave function'
plot 'data/p4_diag.dat'
set term pngcairo
set output "picture/p4_diag.png"
replot
set output
