#!/bin/sh
make -C ../../ mcnde
cp ../../fortran/mcnde .

cat <<EOF > input.nml
&MCPARM
nstochloc = 30,
stochtype = 1,
alpha_s = .1
/

&moistdyn

nonlin_tld = 1.0 !lmd_tld is the sensative param

/
&DATA
TEND    =   800.0
TENERGY =  6.0,
ASST    =  0.200000000000000
toggle_nonlinear = .false.
stochastic_cmt = .false.
dx= 40.0
/
EOF

if [[ -f real.bin ]]
then
    echo "Output file 'real.bin' already exits. Delete it to run model"
    exit -1
fi

cat run.sh > output.log
rm -f column.bin real.bin int.bin OUTPUT snap_shots data.nc && echo 2 | ./mcnde  > output.log 2>&1
./bin2nc.py nc real.bin out.nc
