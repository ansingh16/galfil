f90comp=gfortran
objects = SLICING.o COMB.o
libs= -lblas -llapack
cython=cythonize -a -i





execslice: $(objects)
	$(f90comp) -o execslice $(objects) $(libs)


SLICING.o: SLICING.f90 
	$(f90comp) -c SLICING.f90 $(libs)    

	
COMB.o: COMB.f90
	$(f90comp) -c COMB.f90  $(libs)

%.o: %.f90
	$(f90comp) -c $(libs) $<


clean: 
	rm $(objects)
	
	
	
