set terminal wxt size 1200,900




set palette rgbformulae 22,13,-31
unset key


set multiplot layout 2,2 rowsfirst scale 1,1 title 'Re = 1000 (t = 130 s, uref = 0.1 m/s, 129 x 129)'
# --- GRAPH a


set title 'U'


set xrange [0:1]
set yrange [0:1]
set size ratio 1


set xtics axis
set xtics format '%.1f'
set xlabel ('x [m]')




set ytics axis
set ytics format '%.1f'
set ylabel ('y [m]')


set cblabel ('u [m/s]')




plot 'data/u_129_129_1000_130.txt' with image






# --- GRAPH b


set title 'v'


set xrange [0:1]
set yrange [0:1]
set size ratio 1


set xtics axis
set xtics format '%.1f'
set xlabel ('x [m]')


set ytics axis
set ytics format '%.1f'
set ylabel ('y [m]')


set cblabel ('v [m/s]')






plot 'data/v_129_129_1000_130.txt' with image


# --- GRAPH c




set title 'velocity field'


set xrange [0:1]
set yrange [0:1]
set size ratio 1


set xtics axis
set xtics format '%.1f'
set xlabel ('x [m]')


set ytics axis
set ytics format '%.1f'
set ylabel ('y [m]')


set cblabel ('vel [m/s]')






plot 'data/vel_129_129_1000_130.txt' with image, 'data/vel_129_129_1000_130.txt' using 1:2:($4/(20*sqrt(($4)**2 + ($5)**2))):($5/(20*sqrt(($4)**2 + ($5)**2))) every 5:5 with vectors lc -1 filled






# --- GRAPH d


set title 'pressure'


set xrange [0:1]
set yrange [0:1]
set size ratio 1


set xtics axis
set xtics format '%.1f'
set xlabel ('x [m]')


set ytics axis
set ytics format '%.1f'
set ylabel ('y [m]')


set cblabel ('p [Pa]')






plot 'data/p_129_129_1000_130.txt' with image




unset multiplot
pause 50
