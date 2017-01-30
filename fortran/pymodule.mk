FC= gfortran
FFLAGS = -O3

all: mc.pyf multicloud_mod.o 
	FC=$(FC) f2py -c mc.pyf mc.f90 multicloud_mod.o util.o param_mod.o duni.f ext/minpack/*.f
mc.pyf: mc.f90
	FC=$(FC) f2py -h mc.pyf -m mc mc.f90 --overwrite-signature
.PHONY: mc.pyf

multicloud_mod.o : multicloud_mod.f90 util.o param_mod.o
	$(FC) $(FFLAGS) -fPIC -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -fPIC -c $^

#mcf90.o: multicloud_mod.o param_mod.o util.o duni.f
#	$(FC) $(FFLAGS) -c -o $@ multicloud_mod.o util.o param_mod.o duni.f ext/minpack/*.f

clean: 
	rm -f *.o *.mod
