set terminal wxt

plot "p_129_129_100.txt" with image
set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format "%.1f"
set xlabel ("$x \\ (\\mathrm{m})$")

set ytics axis
set ytics format "%.1f"
set ylabel ("$y \\ (\\mathrm{m})$")

set cblabel ("$u$")

unset key
set palette rgbformulae 22,13,-31
set palette rgb 33,13,10

replot
pause 30
