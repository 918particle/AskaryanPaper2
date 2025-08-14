set terminal pdf enhanced
set output "July7th_plot2.pdf"
unset grid
set xrange [-20:108]
set yrange [-1000:1000]
set xlabel "Time [ns]" font "Courier,18" offset 0,-0.5
set ylabel "Amplitude [mV]" font "Courier,18" offset -2,0
set xtics font "Courier,18"
set ytics font "Courier,18"
set lmargin 10
set rmargin 1
set tmargin 1
set bmargin 4
set key right top box on font "Courier,12" width 1.1 spacing 1.1 at 105,915
set style line 1 lw 1.0 lc rgb "#000000"
set style line 2 lw 2.5 lc rgb "#000000"
set style line 3 lw 5.0 lc rgb "#AAAAAA"
set label 1 "f_0: 0.15 GHz, {/Symbol g}: 0.033 GHz, {/Symbol s}_t: 1.0 ns" at 16,485 font "Courier,14"
set label 2 "R_0 = 1 m ns^{-1}, E_0 = 1 V m^{-1} ns^{-1}" at 16,345 font "Courier,14"
plot "results_July7th_2.dat" using 1:3 w l ls 3 title "mathematical convolution", \
"results_July7th_2.dat" using 1:2 w l ls 1 title "numerical convolution"