set xlabel 'X'
set ylabel 'potential'
plot 'data/potential.dat'
set term pngcairo
set output "picture/potential.png"
replot
set output
