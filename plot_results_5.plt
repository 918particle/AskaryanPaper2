set terminal pdf enhanced
set output "Sept25_plot1.pdf"
unset grid
set xrange [-10:295]
set yrange [-0.05:0.5]
set xlabel "Time [ns]" font "Courier,20" offset 0,-0.5
set ylabel "Normalized Amplitude" font "Courier,20" offset -2,0
set xtics font "Courier,20"
set ytics font "Courier,20"
set lmargin 10
set rmargin 1
set tmargin 1
set bmargin 4
set key right top box on font "Courier,16" width 1.0 spacing 1.25 at 290,0.45
set style line 1 lw 2.5 lc rgb "#000000"
set pointsize 0.5 
plot "example_event.dat" using 1:2 w p pt 6 lc -1 title "Coherently summed waveform (CSW)", "example_event.dat" using 1:3 w l ls 1 title "Analytic envelope ({/symbol r} = 0.951)"
