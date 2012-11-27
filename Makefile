
RCOMPIL = r.compile
RBUILD  = r.build
FLAGS   = -O 2
LFLAGS  = 

OBJ    = cdf2fst.o
OBJ2   = fst_hyb.o

LIBS   = netcdf rmn_013

%.o: %.f90
	$(RCOMPIL) -src $< $(FLAGS)

fst2cdf: $(OBJ)
	$(RBUILD) -obj $^ -arch $(ARCH) -abi $(ABI) $(LFLAGS) -o cdf2fst -libappl "$(LIBS)"

clean:
	rm -f cdf2fst *.o *.mod

