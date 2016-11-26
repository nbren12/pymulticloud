#!/bin/sh
# FC=gfortran CC=gcc CXX=g++ cmake -DCMAKE_BUILD_TYPE=Release .
# make multicloud

make -C fortran -f pymodule.mk
