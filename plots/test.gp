set terminal wxt size 1200,900


set palette rgbformulae 22,13,-31
unset key

set multiplot layout 2,2 rowsfirst scale 1,1
# --- GRAPH a


set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format "%.1f"
set xlabel ("x")

set ytics axis
set ytics format "%.1f"
set ylabel ("y")

set cblabel ("u [m/s]")


plot "u_129_129_100.txt" with image



# --- GRAPH b


set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format "%.1f"
set xlabel ("x")

set ytics axis
set ytics format "%.1f"
set ylabel ("y")

set cblabel ("v [m/s]")



plot "v_129_129_100.txt" with image

# --- GRAPH c


set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format "%.1f"
set xlabel ("x")

set ytics axis
set ytics format "%.1f"
set ylabel ("y")

set cblabel ("vel [m/s]")



# plot "vel_129_129_100.txt" with image, "vel_129_129_100.txt" using 1:2:4:5 every 5:5 with vectors lc -1 filled
plot "vel_129_129_100.txt" with image, "vel_129_129_100.txt" using 1:2:($4/(20*sqrt(($4)**2 + ($5)**2))):($5/(20*sqrt(($4)**2 + ($5)**2))) every 5:5 with vectors lc -1 filled



# --- GRAPH d


set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format "%.1f"
set xlabel ("x")

set ytics axis
set ytics format "%.1f"
set ylabel ("y")

set cblabel ("p [Pa]")



plot "p_129_129_100.txt" with image


unset multiplot

pause 2
