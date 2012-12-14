
RCOMPIL = r.compile
RBUILD  = r.build
FLAGS   = -O 2 -openmp #-optf=-Mipa=pure
LFLAGS  = -openmp #-optf=-Mipa=pure

OBJ    = cdf2fst.o
OBJP   = cdf2fst-ptom.o
OBJE   = cdf2fst-echam.o

LIBS   = netcdf rmn_013

%.o: %.f90
	$(RCOMPIL) -src $< $(FLAGS) 

cdf2fst: $(OBJ)
	$(RBUILD) -obj $^ -arch $(ARCH) -abi $(ABI) $(LFLAGS) -o cdf2fst       -libappl "$(LIBS)"

cdf2fst-ptom: $(OBJP)
	$(RBUILD) -obj $^ -arch $(ARCH) -abi $(ABI) $(LFLAGS) -o cdf2fst-ptom  -libappl "$(LIBS)"
	
cdf2fst-echam: $(OBJE)
	$(RBUILD) -obj $^ -arch $(ARCH) -abi $(ABI) $(LFLAGS) -o cdf2fst-echam -libappl "$(LIBS)"
	
clean:
	rm -f cdf2fst fst2cdf-ptom fst2cdf-echam *.oo *.o *.mod

