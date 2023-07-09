set terminal wxt size 1200,900


set palette rgbformulae 22,13,-31
unset key


set title 'Re = 1000 (t = 120 s, uref = 0.1 m/s, 129 x 129)'


# --- GRAPH c
set title '\textbf{Velocity} (Re = 1000, $t$ = 120 s)'

set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format '%.1f'
set xlabel ("$x \\ [\\mathrm{m}]$")

set ytics axis
set ytics format '%.1f'
set ylabel ("$y \\ [\\mathrm{m}]$")

set cblabel ("$\\| \\mathbf{v} \\| \\ [\\mathrm{m / s}]$")


plot 'data/vel_129_129_1000_120.txt' with image, 'data/vel_129_129_1000_120.txt' using 1:2:($4/(20*sqrt(($4)**2 + ($5)**2))):($5/(20*sqrt(($4)**2 + ($5)**2))) every 10:10 with vectors lc -1 filled

replot

set terminal epslatex standalone color colortext size 12cm,12cm
set rmargin 2
set output "sec_120.tex"
replot
set output

pause 2

