FC= gfortran
FFLAGS = -O3

all: mc.pyf mcf90.o
	f2py -c mc.pyf mc.f90 mcf90.o
mc.pyf: mc.f90
	f2py -h mc.pyf -m mc mc.f90 --overwrite-signature
.PHONY: mc.pyf

mcf90.o: multicloud_mod.f90 param_mod.f90 duni.f
	$(FC) $(FFLAGS) multicloud_mod.f90 param_mod.f90 duni.f ext/minpack/*.f
