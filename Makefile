
RCOMPIL = r.compile
RBUILD  = r.build
FLAGS   = -O 2 -openmp -optf=-Mipa=pure
LFLAGS  = -openmp -optf=-Mipa=pure

OBJ    = cdf2fst.o
OBJ2   = fst_hyb.o

LIBS   = netcdf rmn_013

%.o: %.f90
	$(RCOMPIL) -src $< $(FLAGS) 

fst2cdf: $(OBJ)
	$(RBUILD) -obj $^ -arch $(ARCH) -abi $(ABI) $(LFLAGS) -o cdf2fst -libappl "$(LIBS)"

clean:
	rm -f cdf2fst *.o *.mod

