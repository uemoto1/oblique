#FC=gfortran
FC=ifort
FFLAG=-axMIC-AVX512 -qopenmp -O3

.PHONY: all clean

all: a.out

clean:
	rm a.out *.o *.mod

a.out:  main.o input_parameter.o math_constants.o phys_constants.o em_field.o
	$(FC) $(FFLAG) -o $@ $^

input_parameter.f90: input_parameter.py
	python3 $< > $@

input_parameter.o: input_parameter.f90
	$(FC) $(FFLAG) -c -o $@ $<

math_constants.o: math/math_constants.f90
	$(FC) $(FFLAG) -c -o $@ $<

phys_constants.o: math/phys_constants.f90
	$(FC) $(FFLAG) -c -o $@ $<

em_field.o: rt/em_field.f90 math_constants.o
	$(FC) $(FFLAG) -c -o $@ $<

main.o: main.f90 input_parameter.o em_field.o
	$(FC) $(FFLAG) -c -o $@ $<

