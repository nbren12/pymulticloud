#!/bin/sh
make -C ../../ mcnde
cp ../../fortran/mcnde .

cat <<EOF > input.nml
&MCPARM
nstochloc = 30,
stochtype = 1,
/

&moistdyn

nonlin_tld = 1.0 !lmd_tld is the sensative param

/
&DATA
TEND    =   800.0
TENERGY =  6.0,
ASST    =  0.500000000000000     
toggle_nonlinear = .false.
stochastic_cmt = .false.
/
EOF

rm real.bin
if [[ -f real.bin ]] 
then
    echo "Output file 'real.bin' already exits. Delete it to run model"
    exit -1
fi

rm -f column.bin real.bin int.bin OUTPUT snap_shots data.nc && echo 2 | ./mcnde
./bin2nc.py nc real.bin out.nc
ncview out.nc
