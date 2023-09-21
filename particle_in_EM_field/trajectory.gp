reset
set title "Trajectory of 1000 Particles in E=(0,0,1/2), B=(y/L, x/L, 0)"
set xlabel "x (cs)"
set ylabel "y (cs)"
set zlabel "z (cs)"

set key off #disable legend

#set range of the axis
set xrange [-1000:1000]
set yrange [-1000:1000]
set zrange [-200:10]

#plot each trajectory

splot "particle_1000.dat" i 0 u 3:4:5 w line
replot for [i=1:999] "particle_1000.dat" index i u 3:4:5 w line
