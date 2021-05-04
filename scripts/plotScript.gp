#!/bin/bash

##################################### latestTimeDir ######################################

# ( cd postProcessing/collapsedFields; ls | grep '[0-9]\+$' | sort -n | tail -n1 | while read n;do cp -r "$n" latesTimeDir;done; )

###################################### I/O ###############################################
		
system "mkdir plots"
set terminal postscript eps enhanced

dns 	= "/home/p285464/scripts/DNS"
les 	= "postProcessing/collapsedFields/latesTimeDir"
PLOTS 	= "plots"

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
stats les."/yPlusMean.xy" every ::row1::row1 using col2 nooutput; yPlus	= STATS_max
stats les."/yPlusMean.xy" every ::row2::row2 using col1 nooutput; y	 	= STATS_max

Ub		= 0.1335
nu		= 2e-5
ut		= yPlus*nu/y 		# LES
ut2		= ut*ut
ut_nu 	= ut/nu
UT		= 0.0079			# DNS @ Re_t = 395
UT2		= UT*UT
UT_nu	= UT/nu
print "yPlus = "; print yPlus
print "ut = "; print ut

offset	= 1
ymax	= 6*offset

######################################## Plot Umean #############################################

LES = les."/UMean_X.xy"
DNS = dns."/DNS.txt"

set xlabel 'y/{/Symbol \144}'
set ylabel '<u> / u_{/Symbol t}'
set xrange [0:1]
set yrange [0:30]
plot  DNS u 1:3 			w l ls 10 ti "DNS"
set key autotitle columnhead
replot LES u 1:($2/ut) 		w lp ls 3 ti "LES M1"
set output PLOTS."/1UMean_y.eps"
replot
unset output
set key noautotitle

set xlabel 'y^{/Symbol +}'
set ylabel '<u> / u_{/Symbol t}'
set xrange [1e-1:400]
set yrange [0:30]
set logscale x
plot DNS u 2:3 					w l ls 10 ti "DNS"
set key autotitle columnhead
replot LES u ($1*ut_nu):($2/ut)	w lp ls 3 ti "LES M1"
set output PLOTS."/1UMean_yPl.eps"
replot
unset output
unset logscale x
set key noautotitle

######################################## Plot u_rms noModEff #############################################

DNS = dns."/RMean.txt"

set xlabel 'y/{/Symbol \144}'
set ylabel '( v_{rms} , w_{rms} + 2u_{/Symbol t} , u_{rms} + 4u_{/Symbol t} ) / u_{/Symbol t}'
set xrange [0:1]
set yrange [0:ymax]
plot     DNS u 1:(sqrt(abs($4))) 		w l ls 10 ti "DNS" # vrms
replot   DNS u 1:(sqrt(abs($5)) + 1*offset) 	w l ls 10 ti "" # wrms + ut
replot   DNS u 1:(sqrt(abs($3)) + 2*offset) 	w l ls 10 ti "" # urms + 2ut
set key autotitle columnhead
LESy = les."/UPrime2Mean_YY.xy"
replot LESy u 1:(sqrt(abs($2)/ut2)) 		w lp ls 3 ti "LES"
LESz = les."/UPrime2Mean_ZZ.xy"
replot LESz u 1:(sqrt(abs($2)/ut2) + 1*offset) 	w lp ls 3 ti ""
LESx = les."/UPrime2Mean_XX.xy"
replot LESx u 1:(sqrt(abs($2)/ut2) + 2*offset) 	w lp ls 3 ti ""
set output PLOTS."/2UPrime2Mean_vs_y.eps"
replot
unset output
set key noautotitle


set xlabel 'y/{/Symbol \144}'
set ylabel '( v_{rms} , w_{rms} + 2u_{/Symbol t} , u_{rms} + 4u_{/Symbol t} ) / u_{/Symbol t}'
set xrange [1e-1:400]
set yrange [0:ymax]
set logscale x
plot     DNS u 2:(sqrt(abs($4))) 		w l ls 10 ti "DNS" # vrms
replot   DNS u 2:(sqrt(abs($5)) + 1*offset) 	w l ls 10 ti "" # wrms + ut
replot   DNS u 2:(sqrt(abs($3)) + 2*offset) 	w l ls 10 ti "" # urms + 2ut
set key autotitle columnhead
LESy = les."/UPrime2Mean_YY.xy"
replot LESy u ($1*ut_nu):(sqrt(abs($2)/ut2)) 		w lp ls 3 ti "LES"
LESz = les."/UPrime2Mean_ZZ.xy"
replot LESz u ($1*ut_nu):(sqrt(abs($2)/ut2) + 1*offset) 	w lp ls 3 ti ""
LESx = les."/UPrime2Mean_XX.xy"
replot LESx u ($1*ut_nu):(sqrt(abs($2)/ut2) + 2*offset) 	w lp ls 3 ti ""
set output PLOTS."/2UPrime2Mean_vs_yPlus.eps"
replot
unset logscale x
unset output
set key noautotitle


######################################## Plot uv noModEff #############################################

DNS = dns."/RMean.txt"
LES = les."/UPrime2Mean_XY.xy"


set xlabel 'y/{/Symbol \144}'
set ylabel '- < uv > / u_{/Symbol t}^2'
set xrange [0:1]
set yrange [0:1]
plot   DNS u 1:(-$6) 		w l ls 10 ti "DNS"
set key autotitle columnhead
replot LES u 1:((-$2)/ut2) 	w lp ls 3 ti "LES"
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
replot LES u ($1*ut_nu):((-$2)/ut2)	w lp ls 3 ti "LES"
set output PLOTS."/3uvUPrime2Mean_vs_yPlus.eps"
replot
unset logscale x
unset output
set key noautotitle

######################################## Plot u_rms withModEff #############################################

DNS = dns."/DNS.txt"

set xlabel 'y/{/Symbol \144}'
set ylabel '( v_{rms} , w_{rms} + 2u_{/Symbol t} , u_{rms} + 4u_{/Symbol t} ) / u_{/Symbol t}'
set xrange [0:1]
set yrange [0:ymax]
plot     DNS u 1:(sqrt(abs($6))) 		w l ls 10 ti "DNS" # vrms
replot   DNS u 1:(sqrt(abs($7)) + 1*offset) 	w l ls 10 ti "" # wrms + ut
replot   DNS u 1:(sqrt(abs($5)) + 2*offset) 	w l ls 10 ti "" # urms + 2ut
set key autotitle columnhead
LESy = les."/RdevPlusModEffMean_YY.xy"
replot LESy u 1:(sqrt(abs($2)/ut2)) 		w lp ls 3 ti "LES"
LESz = les."/RdevPlusModEffMean_ZZ.xy"
replot LESz u 1:(sqrt(abs($2)/ut2) + 1*offset) 	w lp ls 3 ti ""
LESx = les."/RdevPlusModEffMean_XX.xy"
replot LESx u 1:(sqrt(abs($2)/ut2) + 2*offset) 	w lp ls 3 ti ""
set output PLOTS."/4RdevPlusModEffMean_vs_y.eps"
replot
unset output
set key noautotitle


set xlabel 'y/{/Symbol \144}'
set ylabel '( v_{rms} , w_{rms} + 2u_{/Symbol t} , u_{rms} + 4u_{/Symbol t} ) / u_{/Symbol t}'
set xrange [1e-1:400]
set yrange [0:ymax]
set logscale x
plot     DNS u 2:(sqrt(abs($6))) 		w l ls 10 ti "DNS" # vrms
replot   DNS u 2:(sqrt(abs($7)) + 1*offset) 	w l ls 10 ti "" # wrms + ut
replot   DNS u 2:(sqrt(abs($5)) + 2*offset) 	w l ls 10 ti "" # urms + 2ut
set key autotitle columnhead
LESy = les."/RdevPlusModEffMean_YY.xy"
replot LESy u ($1*ut_nu):(sqrt(abs($2)/ut2)) 		w lp ls 3 ti "LES"
LESz = les."/RdevPlusModEffMean_ZZ.xy"
replot LESz u ($1*ut_nu):(sqrt(abs($2)/ut2) + 1*offset) 	w lp ls 3 ti ""
LESx = les."/RdevPlusModEffMean_XX.xy"
replot LESx u ($1*ut_nu):(sqrt(abs($2)/ut2) + 2*offset) 	w lp ls 3 ti ""
set output PLOTS."/4RdevPlusModEffMean_vs_yPlus.eps"
replot
unset logscale x
unset output
set key noautotitle


######################################## Plot uv withModEff #############################################

DNS = dns."/DNS.txt"
LES = les."/RdevPlusModEffMean_XY.xy"


set xlabel 'y/{/Symbol \144}'
set ylabel '- < uv > / u_{/Symbol t}^2'
set xrange [0:1]
set yrange [0:1]
plot   DNS u 1:(-$8) 		w l ls 10 ti "DNS"
set key autotitle columnhead
replot LES u 1:((-$2)/ut2) 	w lp ls 3 ti "LES"
set output PLOTS."/5uvRdevPlusModEffMean_vs_y.eps"
replot
unset output
set key noautotitle

set xlabel 'y^{/Symbol +}'
set ylabel '- < uv > / u_{/Symbol t}^2'
set xrange [1e-1:400]
set yrange [0:1]
set logscale x
plot DNS u 2:(-$8)				w l ls 10 ti "DNS"
set key autotitle columnhead
replot LES u ($1*ut_nu):((-$2)/ut2)	w lp ls 3 ti "LES"
set output PLOTS."/5uvRdevPlusModEffMean_vs_yPlus.eps"
replot
unset logscale x
unset output
set key noautotitle

###################################### Endline ###############################################
