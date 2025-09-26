set terminal pdf enhanced
set output "Sept25_plot2.pdf"
unset grid
set xrange [-10:266]
set yrange [-0.05:0.5]
set xlabel "Time [ns]" font "Courier,18" offset 0,-0.5
set ylabel "Normalized Amplitude" font "Courier,18" offset -2,0
set xtics font "Courier,18"
set ytics font "Courier,18"
set lmargin 10
set rmargin 1
set tmargin 1
set bmargin 4
set key right top box on font "Courier,12" width 1.1 spacing 1.25 at 250,0.45
set style line 1 lw 2.5 lc rgb "#000000"
set pointsize 0.5 
plot "example_event_bad.dat" using 1:2 w lp pt 6 lc -1 title "Coherently summed waveform (CSW)", "example_event_bad.dat" using 1:3 w l ls 1 title "Analytic envelope ({/symbol r} = 0.773)"