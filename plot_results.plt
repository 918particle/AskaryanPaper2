set terminal pdf enhanced
set output "July3rd_plot1.pdf"
unset grid
set xrange [-20:108]
set yrange [-300:300]
set xlabel "Time [ns]" font "Courier,18" offset 0,-0.5
set ylabel "Amplitude [mV]" font "Courier,18" offset -2,0
set xtics font "Courier,18"
set ytics font "Courier,18"
set lmargin 10
set rmargin 1
set tmargin 1
set bmargin 4
set key right top box on font "Courier,12" width 1.1 spacing 1.1 at 105,275
set style line 1 lw 1 lc rgb "#000000"
set style line 2 lw 2 lc rgb "#444444"
set style line 3 lw 4 lc rgb "#AAAAAA"
set multiplot
set origin 0,0
set size 1,1
set label 1 "f_0: 0.30 GHz, {/Symbol g}: 0.05 GHz, {/Symbol s}_t: 0.5 ns" at 16,145 font "Courier,10.5"
set label 2 "R_0 = 1 m ns^{-1}, E_0 = 1 V m^{-1} ns^{-1}" at 16,115 font "Courier,10.5"
plot "results_July3rd.dat" using 1:2 w l ls 1 title "mathematical convolution", \
"results_July3rd.dat" using 1:3 w l ls 3 title "envelope of mathematical convolution", \
"results_July3rd.dat" using 1:4 w l ls 2 title "mathematical envelope"
set origin 0.4,0.1
set size 0.55,0.45
set xrange [-5:20]
set yrange [-300:300]
unset key
unset xlabel
unset ylabel
unset label 1
unset label 2
set xtics font "Courier,14"
set ytics font "Courier,14"
plot "results_July3rd.dat" using 1:2 w l ls 1 title "mathematical convolution", \
"results_July3rd.dat" using 1:3 w l ls 3 title "envelope of mathematical convolution", \
"results_July3rd.dat" using 1:4 w l ls 2 title "mathematical envelope"