# ################################################################## #
# Usage: gnuplot -e "p='01'" plot_block_comparison.gnu
# ################################################################## #
set terminal postscript eps enhanced color
set output 'block_compare_'.p.'.eps'

set grid

set key left top

set title "Average time to access random position (p = ".p[1:1].".".p[2:].")"

set style line 1 linetype 1 linecolor rgb "black" pt 6
set style line 2 linetype 2 linecolor rgb "black" pt 3

set xlabel 'Block\_size [bits]'
set ylabel 'Access time [CPU cycles]'

plot 'block_access_'.p.'.data' u 1:4  title 'RRR'  with linespoint linestyle 1, \
     'block_access_'.p.'.data' u 1:6  title 'R3D3' with linespoint linestyle 2
