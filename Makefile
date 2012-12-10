
RCOMPIL = r.compile
RBUILD  = r.build
FLAGS   = -O 2 -openmp -optf=-Mipa=pure
LFLAGS  = -openmp -optf=-Mipa=pure

OBJ    = cdf2fst.o
OBJP   = cdf2fst-ptom.o

LIBS   = netcdf rmn_013

%.o: %.f90
	$(RCOMPIL) -src $< $(FLAGS) 

fst2cdf: $(OBJ)
	$(RBUILD) -obj $^ -arch $(ARCH) -abi $(ABI) $(LFLAGS) -o cdf2fst -libappl "$(LIBS)"

fst2cdf-ptom: $(OBJP)
	$(RBUILD) -obj $^ -arch $(ARCH) -abi $(ABI) $(LFLAGS) -o cdf2fst-ptom -libappl "$(LIBS)"
	
clean:
	rm -f cdf2fst *.o *.mod

