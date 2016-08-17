#!/bin/sh
make

cat <<EOF > input.nml
&DATA
TEND    =   100.0
TENERGY =  1.0,
ASST    =  0.000000000000000     
toggle_nonlinear = .true.
stochastic_cmt = .true.
/
EOF

rm -f column.bin real.bin int.bin OUTPUT snap_shots data.nc && echo 2 | ./column
python bin2nc.py column column.bin
