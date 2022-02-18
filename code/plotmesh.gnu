# Type of terminal (wxt, x11, qt, pngcairo, ...) and size
set terminal qt size 700,700

# Line styles
set style line 1 lc rgb 'black' lw 2        # Domain border
set style line 2 lc rgb 'blue' pt 7 ps 2    # Nodes
set style line 3 lc rgb 'red' lw 1 dt 2     # Control volume

# No legend
set nokey

# Plot the grid
plot "data/domain_data.dat" using 1:2:($3-$1):($4-$2) with vectors nohead ls 1, \
    "data/facex_data.dat" using 1:2:($3-$1):($4-$2) with vectors nohead ls 3, \
    "data/facey_data.dat" using 1:2:($3-$1):($4-$2) with vectors nohead ls 3, \
    "data/node_data.dat" with points ls 2, \

# Same scale for x and y axis
set size ratio -1

# Axis labels
set xlabel "x [m]"
set ylabel "y [m]"

# Prevent window from closing
pause -1
