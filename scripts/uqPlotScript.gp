#!/bin/bash

###################################### I/O ###############################################
		
system "mkdir plots"
set terminal postscript eps enhanced
set encoding iso_8859_1 

dns 	= "/home/p285464/scripts/DNS"
les 	= "postProcessing/collapsedFields/latesTimeDir"
YPlus	= "postProcessing/patchExpression_yPlus/0/bottomWall"
PLOTS 	= "plots"
LESy 	= les."/pMean.xy"

###################################### Style ###############################################

set style line 1  lc rgb '#00ff00' lw 2 lt 1 pt 1 #green
set style line 2  lc rgb '#0000ff' lw 2 lt 1 pt 1 #blue
set style line 3  lc rgb '#ff0000' lw 2 lt 1 pt 7 #red
set style line 10 lc rgb 'black'   lw 2 lt 1 pt 1 #green

set grid
set tics out nomirror

######################################## Initialize #############################################

#### Constants ####
SMALL = 1e-12

#### Find u_tau ####

row1 = 0; row2 = 1; col1 = 1; col2 = 2
stats YPlus using col2 nooutput; 					yPlus	= STATS_mean
stats LESy every ::row2::row2 using col1 nooutput; 	y_wall	= STATS_max

Ub		= 0.1335
nu		= 2e-5
#ut		= 1.8*nu/y_wall 
ut		= yPlus*nu/y_wall 		# LES
ut2		= ut*ut
ut_nu 	= ut/nu
UT		= 0.0079			# DNS @ Re_t = 395
UT2		= UT*UT
UT_nu	= UT/nu

N		= 1 # times sigma
print "yPlus = "; print yPlus
print "ut = "; print ut

######################################## Plot pMean #############################################

LES 	= les."/pMean.xy"

set xlabel 'y/{/Symbol \144}'
set ylabel '< p >'
set xrange [0:1]
set yrange [-0.02:0.02]
set key autotitle columnhead
plot \
LES u 1:($2/ut) 			w lp ls 3 ti "LES", \
LES u 1:(($2 + N*$3)/ut) 	w l ls 1 ti "( < u > \261 ".N."{/Symbol s}_{< u >} ) / u_{/Symbol t}", \
LES u 1:(($2 - N*$3)/ut) 	w l ls 1 ti ""
set output PLOTS."/0pMean_vs_y.eps"
replot
unset output
set key noautotitle

######################################## Plot UMean #############################################

LES 	= les."/UMean_X.xy"
DNS 	= dns."/DNS.txt"

set xlabel 'y/{/Symbol \144}'
set ylabel '< u > / u_{/Symbol t}'
set xrange [0:1]
set yrange [0:30]
plot  DNS u 1:3 		w l ls 10 ti "DNS"
set key autotitle columnhead
replot \
LES u 1:($2/ut) 			w lp ls 3 ti "LES", \
LES u 1:(($2 + N*$3)/ut) 	w l ls 1 ti "( < u > \261 ".N."{/Symbol s}_{< u >} ) / u_{/Symbol t}", \
LES u 1:(($2 - N*$3)/ut) 	w l ls 1 ti ""
set output PLOTS."/1UMean_vs_y.eps"
replot
unset output
set key noautotitle

set xlabel 'y^{/Symbol +}'
set ylabel '< u > / u_{/Symbol t}'
set xrange [1e-1:400]
set yrange [0:30]
set logscale x
plot DNS u 2:3 					w l ls 10 ti "DNS"
set key autotitle columnhead
replot \
LES u ($1*ut_nu):($2/ut)		w lp ls 3 ti "LES", \
LES u ($1*ut_nu):(($2 + N*$3)/ut) w l ls 1 ti "( < u > \261 ".N."{/Symbol s}_{< u >} ) / u_{/Symbol t}", \
LES u ($1*ut_nu):(($2 - N*$3)/ut) w l ls 1 ti ""
set output PLOTS."/1UMean_vs_yPl.eps"
replot
unset output
unset logscale x
set key noautotitle

######################################## Plot u_rms noModEff #############################################

DNS = dns."/RMean.txt"

set xlabel 'y/{/Symbol \144}'
set ylabel '( v_{rms} , w_{rms} + 2u_{/Symbol t} , u_{rms} + 4u_{/Symbol t} ) / u_{/Symbol t}'
set xrange [0:1]
set yrange [0:10]
plot \
DNS u 1:(sqrt(abs($4))) 		w l ls 10 ti "DNS", \
DNS u 1:(sqrt(abs($5)) + 2) 	w l ls 10 ti "", \
DNS u 1:(sqrt(abs($3)) + 4) 	w l ls 10 ti "" 	 	
set key autotitle columnhead
LESy = les."/UPrime2Mean_YY.xy"
replot \
LESy u 1:(sqrt(abs($2)/ut2)) 		w lp ls 3 ti "LES", \
LESy u 1:((sqrt(abs($2)) + N*sqrt(abs($3)))/ut) w l ls 1 ti "( u_{rms} \261 ".N."{/Symbol s}_{u_{rms}} ) / u_{/Symbol t}", \
LESy u 1:((sqrt(abs($2)) - N*sqrt(abs($3)))/ut) w l ls 1 ti ""
LESz = les."/UPrime2Mean_ZZ.xy"
replot \
LESz u 1:(sqrt(abs($2)/ut2) + 2) 	w lp ls 3 ti "", \
LESz u 1:((sqrt(abs($2)) + N*sqrt(abs($3)))/ut + 2) w l ls 1 ti "", \
LESz u 1:((sqrt(abs($2)) - N*sqrt(abs($3)))/ut + 2) w l ls 1 ti ""
LESx = les."/UPrime2Mean_XX.xy"
replot \
LESx u 1:(sqrt(abs($2)/ut2) + 4) 	w lp ls 3 ti "", \
LESx u 1:((sqrt(abs($2)) + N*sqrt(abs($3)))/ut + 4) w l ls 1 ti "", \
LESx u 1:((sqrt(abs($2)) - N*sqrt(abs($3)))/ut + 4) w l ls 1 ti ""
set output PLOTS."/2UPrime2Mean_vs_y.eps"
replot
unset output
set key noautotitle


set xlabel 'y/{/Symbol \144}'
set ylabel '( v_{rms} , w_{rms} + 2u_{/Symbol t} , u_{rms} + 4u_{/Symbol t} ) / u_{/Symbol t}'
set xrange [1e-1:400]
set yrange [0:10]
set logscale x
plot \
DNS u 2:(sqrt(abs($4))) 		w l ls 10 ti "DNS", \
DNS u 2:(sqrt(abs($5)) + 2) 	w l ls 10 ti "", \
DNS u 2:(sqrt(abs($3)) + 4) 	w l ls 10 ti ""
set key autotitle columnhead
LESy = les."/UPrime2Mean_YY.xy"
replot \
LESy u ($1*ut_nu):(sqrt(abs($2)/ut2)) 		w lp ls 3 ti "LES", \
LESy u ($1*ut_nu):((sqrt(abs($2)) + N*sqrt(abs($3)))/ut) w l ls 1 ti "( u_{rms} \261 ".N."{/Symbol s}_{u_{rms}} ) / u_{/Symbol t}",\
LESy u ($1*ut_nu):((sqrt(abs($2)) - N*sqrt(abs($3)))/ut) w l ls 1 ti ""
LESz = les."/UPrime2Mean_ZZ.xy"
replot \
LESz u ($1*ut_nu):(sqrt(abs($2)/ut2) + 2) 	w lp ls 3 ti "", \
LESz u ($1*ut_nu):((sqrt(abs($2)) + N*sqrt(abs($3)))/ut + 2) w l ls 1 ti "", \
LESz u ($1*ut_nu):((sqrt(abs($2)) - N*sqrt(abs($3)))/ut + 2) w l ls 1 ti ""
LESx = les."/UPrime2Mean_XX.xy"
replot \
LESx u ($1*ut_nu):(sqrt(abs($2)/ut2) + 4) 	w lp ls 3 ti "", \
LESx u ($1*ut_nu):((sqrt(abs($2)) + N*sqrt(abs($3)))/ut + 4) w l ls 1 ti "", \
LESx u ($1*ut_nu):((sqrt(abs($2)) - N*sqrt(abs($3)))/ut + 4) w l ls 1 ti ""
set output PLOTS."/2UPrime2Mean_vs_yPlus.eps"
replot
unset logscale x
unset output
set key noautotitle

######################################## Plot uv noModEff #############################################

DNS = dns."/RMean.txt"
LES = les."/UPrime2Mean_XY.xy"
#N = 1

set xlabel 'y/{/Symbol \144}'
set ylabel '- < uv > / u_{/Symbol t}^2'
set xrange [0:1]
set yrange [0:1]
plot   DNS u 1:(-$6) 		  w l ls 10 ti "DNS"
set key autotitle columnhead
replot \
LES u 1:((-$2)/ut2) 		  w lp ls 3 ti "LES", \
LES u 1:(((-$2) + N*(-$3))/ut2) w l ls 1 ti "- ( < uv > \261 ".N."{/Symbol s}_{< uv >} ) / u_{/Symbol t}", \
LES u 1:(((-$2) - N*(-$3))/ut2) w l ls 1 ti ""
set output PLOTS."/3uvUPrime2Mean_vs_y.eps"
replot
unset output
set key noautotitle

set xlabel 'y^{/Symbol +}'
set ylabel '- < uv > / u_{/Symbol t}^2'
set xrange [1e-1:400]
set yrange [0:1]
set logscale x
plot DNS u 2:(-$6)				w l ls 10 ti "DNS"
set key autotitle columnhead
replot \
LES u ($1*ut_nu):((-$2)/ut2)	w lp ls 3 ti "LES", \
LES u ($1*ut_nu):(((-$2) + N*(-$3))/ut2) w l ls 1 ti "- ( < uv > \261 ".N."{/Symbol s}_{< uv >} ) / u_{/Symbol t}", \
LES u ($1*ut_nu):(((-$2) - N*(-$3))/ut2) w l ls 1 ti ""
set output PLOTS."/3uvUPrime2Mean_vs_yPlus.eps"
replot
unset logscale x
unset output
set key noautotitle

###################################### Endline ###############################################
