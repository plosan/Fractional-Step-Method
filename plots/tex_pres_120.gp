set terminal wxt size 1200,900


set palette rgbformulae 22,13,-31
unset key


set title 'Re = 1000 (t = 120 s, uref = 0.1 m/s, 129 x 129)'

# --- GRAPH d
set title '\textbf{Pressure} (Re = 1000, $t$ = 120 s)'

set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format '%.1f'
set xlabel ("Position $x \\ [\\mathrm{m}]$")

set ytics axis
set ytics format '%.1f'
set ylabel ("Position $y \\ [\\mathrm{m}]$")

set cbrange [-0.002:0.01]
set cbtics format '%.3f'
set cblabel ("Pressure - $10^5$ $[\\mathrm{Pa}]$")

plot '../output_data/p_129_129_1000_120.txt' using 1:2:($3 - 100000) with image

replot

set terminal epslatex standalone color colortext size 12cm,12cm
set rmargin 2
set output "pres_120.tex"
replot
set output

pause 2