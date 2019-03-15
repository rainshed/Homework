set xlabel 'X'
set ylabel 'Wave function'
plot 'data/p3n_wave_function.dat'
set term pngcairo
set output "picture/p3_numerov.png"
replot
set output
