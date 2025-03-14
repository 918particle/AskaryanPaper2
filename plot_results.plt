set terminal pdf enhanced
set output "March12_plot1.pdf"
unset grid
set xrange [-20:108]
set yrange [-100:100]
set xlabel "Time [ns]" font "Courier,18" offset 0,-0.5
set ylabel "Amplitude [mV]" font "Courier,18" offset -2,0
set xtics font "Courier,18"
set ytics font "Courier,18"
set lmargin 10
set rmargin 1
set tmargin 1
set bmargin 4
set key right top box on font "Courier,12" width 1.1 spacing 1.1 at 105,90
set style line 1 lw 0.75 lc rgb "#000000"
set style line 2 lw 1.5 lc rgb "#FFCC00"
set style line 3 lw 3 lc rgb "#5500AA"
set multiplot
set origin 0,0
set size 1,1
set label 1 "f_0: 0.55 GHz, {/Symbol g}: 0.4 GHz, {/Symbol s}_t: 0.25 ns" at 40,50 font "Courier,10.5"
set label 2 "R_0 = 1 m ns^{-1}, E_0 = 1 V m^{-1} ns^{-1}" at 40,37 font "Courier,10.5"
plot "results_March12.dat" using 1:4 w l ls 3 title "mathematical envelope", \
"results_March12.dat" using 1:3 w l ls 2 title "envelope of mathematical convolution", \
"results_March12.dat" using 1:2 w l ls 1 title "mathematical convolution"
set origin 0.4,0.1
set size 0.55,0.45
set xrange [8:31]
set yrange [-100:100]
unset key
unset xlabel
unset ylabel
unset label 1
unset label 2
set xtics font "Courier,14"
set ytics font "Courier,14"
plot "results_March12.dat" using 1:4 w l ls 3 title "mathematical envelope", \
"results_March12.dat" using 1:3 w l ls 2 title "envelope of mathematical convolution", \
"results_March12.dat" using 1:2 w l ls 1 title "mathematical convolution"