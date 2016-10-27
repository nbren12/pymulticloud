#!/bin/sh
make -C ../../ mcnde
cp ../../fortran/mcnde .

cat <<EOF > input.nml
&MCPARM
nstochloc = 60,
stochtype = 1,
theta_eb_m_theta_em = 11.0,
theta_ebs_m_theta_eb = 15,
xis= 0.00, !doesn't matter
mu = .00,  ! doesn't matter
alpha3=0.0,  ! doesn't matter
a1 = 0.0, ! a1 vs a2 does not matter much although moisture closure gives fishbone waves
deltac1= 0.0,  ! congestus detrainment matters a little. but not much
alpha2 = 0.0,  ! doesn't matter
alpha3 = 0.0,  ! doesn't matter
alpha4 = 4.0,  ! matters a lot
! tau30=1.0      ! matters a lot
alpha_s=0.25      ! matters a lot

/

&moistdyn

nonlin_tld = 0.0 !lmd_tld is the sensative param

/
&DATA
TEND    =   50.0
TENERGY =  6.0,
ASST    =  0.300000000000000     
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
