set xlabel 'X'
set ylabel 'Wave function'
plot 'data/p4_wave_function.dat'
set term pngcairo
set output "picture/p4_wave_function.png"
replot
set output
