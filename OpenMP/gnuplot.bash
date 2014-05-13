#!/usr/bin/gnuplot

#runs commands for gnuplot
reset
set terminal png
set xlabel "Number of Crystals"
set ylabel "Time (Seconds)"
set title "OMP Scaling"  
set term png
set output "OMP_Scaling.png"

#set   autoscale                        # scale axes automatically

                        # set xtics automatically
#set ytic auto                          # set ytics automatically
 

plot "timing.txt" using 1:2 with linespoints
set logscale x 2
