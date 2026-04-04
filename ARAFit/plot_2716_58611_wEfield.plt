set terminal pdf enhanced
set output "2716_58611_wEfield.pdf"
unset grid
set xrange [75:275]
set yrange [-0.5:0.5]
set xlabel "Time [ns]" font "Courier,18" offset 0,-0.5
set ylabel "Amplitude [norm.]" font "Courier,18" offset -2,0
set xtics 75,50,275 font "Courier,18"
set ytics font "Courier,18"
set lmargin 10
set rmargin 5
set tmargin 1
set bmargin 4
set key right top box on font "Courier,14" width 1.1 spacing 1.1 at 262.5,0.45
set style line 1 lw 2.0 lc rgb "#AAAAAA" 
set style line 2 lw 1.5 lc rgb "#000000" pt 6 ps 0.45
set style line 3 lw 2.5 lc rgb "#000000" pt 6 ps 0.45
set style line 4 lw 2.5 lc rgb "#AAAAAA" pt 6 ps 0.45
set style line 5 lw 0.75 lc rgb "#000000" pt 6 ps 0.45
set multiplot
set origin 0,0
set size 0.5,1
plot "test.dat" using 1:2 w l ls 1 title "CSW", "test.dat" using 1:3 w l ls 2 title "Model"
set origin 0.5,0
set size 0.5,1
set xlabel "Time [ns]" font "Courier,18" offset 0,-0.5
set ylabel "E-field [norm.]" font "Courier,18" offset -2,0
unset label 1
unset label 2
unset label 3
unset label 4
set xrange [100:125]
set yrange [-1.0:1.25]
set xtics 100,10,125 font "Courier,18"
set ytics font "Courier,18"
set key right top box on font "Courier,14" width 1.1 spacing 1.1 at 123.5,1.17
plot "test.dat" using 1:4 w l ls 3 title "{/Symbol s}=0.6 ns", \
"test.dat" using 1:5 w l ls 4 title "{/Symbol s}-0.1 ns", \
"test.dat" using 1:6 w l ls 4 title "{/Symbol s}+0.1 ns", \
"adapted.dat" using ($1+111.5):(-$2) w l ls 5 title "UHECRs"