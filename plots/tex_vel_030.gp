set terminal wxt size 1200,900


set palette rgbformulae 22,13,-31
unset key


set title 'Re = 1000 (t = 30 s, uref = 0.1 m/s, 129 x 129)'


# --- GRAPH c
set title '\textbf{Velocity field} (Re = 1000, $t$ = 30 s)'

set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format '%.2f'
set xlabel ("Position $x \\ [\\mathrm{m}]$")

set ytics axis
set ytics format '%.2f'
set ylabel ("Position $y \\ [\\mathrm{m}]$")

set cblabel ("Velocity $[\\mathrm{m / s}]$")
set cbtics format '%.2f'


plot '../output_data/vel_129_129_1000_30.txt' with image, '../output_data/vel_129_129_1000_30.txt' using 1:2:($4/(20*sqrt(($4)**2 + ($5)**2))):($5/(20*sqrt(($4)**2 + ($5)**2))) every 10:10 with vectors lc -1 filled

replot

set terminal epslatex standalone color colortext size 15cm,12cm
set rmargin 2
set output "vel_030.tex"
replot
set output

pause 2

