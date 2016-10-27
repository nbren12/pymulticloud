#!/bin/sh
make -C ../../ mcnde
cp ../../fortran/mcnde .

cat <<EOF > input.nml
&MCPARM
nstochloc = 60,
stochtype = 1,
taumult = .25,
/
&DATA
TEND    =   400.0
TENERGY =  6.0,
ASST    =  0.000000000000000     
toggle_nonlinear = .false.
stochastic_cmt = .false.
/
EOF

if [[ -f real.bin ]] 
then
    echo "Output file 'real.bin' already exits. Delete it to run model"
    exit -1
fi

rm -f column.bin real.bin int.bin OUTPUT snap_shots data.nc && echo 2 | ./mcnde
python bin2nc.py report real.bin
