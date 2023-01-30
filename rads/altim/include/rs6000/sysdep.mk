F90     = xlf
FFLAGS 	= -u -qmaxmem=-1
RECUR   = -qrecur
CPLUS	= -qcpluscmt
MATHLIB = -L/usr/local/lib -llapack
SHARED_LD       = buildshared
SHARED_LIB_LIBS = -lX11 -lxlf90 -lm -lc
SHARED_EXT	= so

#.F.o mechanism
.SUFFIXES:	.F
.F.o:
	$(CPP) $< > F$*.f
	$(FC) $(FFLAGS) -c F$*.f
	mv F$*.o $@
	rm -f F$*.f
