#!/usr/bin/gnuplot

#runs commands for gnuplot
reset

# Fit data
f(x) = m*x
fit f(x) 'timing_TR1024.txt' using 1:2 via m

set terminal png
set xlabel "Number of Crystals"
set ylabel "Time (Seconds)"
set title "Spectral Scaling"  
set term png
set output "Spectral_Scaling.png"
                      
set ytic auto           # set ytics automatically
set xtic auto		# set xtics automatically
#set logscale y 10
#set logscale x 2
set key left top
plot "timing_TR1024.txt" using 1:2 with linespoints title 'Timing', f(x) title 'f=m*x'


