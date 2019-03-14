set term pngcairo
set xlabel 'X'
set ylabel 'Wave function'
set output "picture/p1_wave_function1.png"
plot 'data/p1_wave_function1.dat'
set output
set output "picture/p1_wave_function2.png"
plot 'data/p1_wave_function2.dat'
set output
set output "picture/p1_wave_function3.png"
plot 'data/p1_wave_function3.dat'
set output
