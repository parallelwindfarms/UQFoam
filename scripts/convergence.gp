#!/bin/bash

###################################### I/O ###############################################
		
system "mkdir plots"
set terminal postscript eps enhanced
set encoding iso_8859_1 

PLOTS = "plots"

###################################### Style ###############################################

set grid
set tics out nomirror
set key autotitle columnhead

######################################## Plot pMean #############################################

set xlabel 't (s)'
set ylabel '< p >'
plot "0/pMean0" u 1:($2) w l ti ""
set output PLOTS."/pMean.eps"
replot
unset output

######################################## Plot UMean #############################################

set xlabel 't (s)'
set ylabel '< u >'
#set yrange [0.1:0.2]
plot "0/UMean0" u 1:($2) " %lf (%lf %lf %lf)" w l
set output PLOTS."/UMeanConv.eps"
replot
unset output

######################################## Plot USigmaMean #############################################

set xlabel 't (s)'
set ylabel '{/Symbol s}_{< u >}'
#set yrange [0:2e-4]
plot "0/UMeanSigma" u 1:($2) " %lf (%lf %lf %lf)" w l
set output PLOTS."/USigmaMeanConv.eps"
replot
unset output

