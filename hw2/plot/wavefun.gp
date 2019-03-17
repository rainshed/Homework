set xlabel 'X'
set ylabel 'Wave function'

set title 'wave funtion the ground state'
plot 'data/p1_wave_function1.dat'
pause 2

set title 'wave funtion of the 1st excited state' 
plot 'data/p1_wave_function2.dat'
pause 2

set title 'wave funtion of the 2st excited state' 
plot 'data/p1_wave_function3.dat'
pause 2

set term pngcairo

set output "picture/p1_wave_function1.png"
set title 'wave funtion the ground state'
plot 'data/p1_wave_function1.dat'
set output

set output "picture/p1_wave_function2.png"
set title 'wave funtion of the 1st excited state' 
plot 'data/p1_wave_function2.dat'
set output

set output "picture/p1_wave_function3.png"
set title 'wave funtion of the 2st excited state' 
plot 'data/p1_wave_function3.dat'
set output
