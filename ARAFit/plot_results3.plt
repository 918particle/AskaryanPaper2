set terminal pdf enhanced
set output "November10th_plot2.pdf"
unset grid
set xrange [120:340]
set yrange [-0.5:0.5]
set xlabel "Time [ns]" font "Courier,18" offset 0,-0.5
set ylabel "Amplitude [norm.]" font "Courier,18" offset -2,0
set xtics 125,25,350 font "Courier,18"
set ytics font "Courier,18"
set lmargin 10
set rmargin 5
set tmargin 1
set bmargin 4
set key right top box on font "Courier,14" width 1.1 spacing 1.1 at 325,0.45
set style line 1 lw 3.0 lc rgb "#AAAAAA"
set style line 2 lw 1.0 lc rgb "#000000"
set label 1 "{/Symbol s}_t: 0.425 ns" at 230,0.3 font "Courier,12"
set label 2 "f_0: 0.195 GHz" at 230,0.25 font "Courier,12"
set label 3 "{/Symbol g}: 0.019 GHz" at 230,0.2 font "Courier,12"
set label 4 "{/Symbol r}: 0.770" at 230,0.15 font "Courier,12"
set multiplot
set origin 0,0
set size 1,1
plot "result_Nov10_2.dat" using 1:2 w l ls 1 title "CSW, ARA CR event", \
"result_Nov10_2.dat" using 1:3 w l ls 2 title "Model"
set origin 0.4,0.1
set size 0.55,0.45
set xrange [0:350]
set yrange [-0.5:0.5]
unset key
unset xlabel
unset ylabel
unset label 1
unset label 2
unset label 3
unset label 4
set xtics 0,100,350  font "Courier,14"
set ytics font "Courier,14"
plot "result_Nov10_2.dat" using 1:2 w l ls 1, "result_Nov10_2.dat" using 1:3 w l ls 2