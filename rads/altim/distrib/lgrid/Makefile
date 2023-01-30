########################################################################
#
#	Makefile for LGRID
#
# 26-Feb-1998 - Remko Scharroo - Multi-platform version
#
########################################################################
#
# This Makefile contains the description for the generation and testing of
# a program to convert the DEOS binary grid format into ASCII.
#
# Run 'make' to generate the executable 'lgrid', run 'make test' to test
# its functioning.
#
# It has been attempted to make this Makefile platform independent, but of
# course it might be necessary to change a few things in the following
# definitions in case the compilation or the test was unsuccessful.
#
# For information please contact Remko Scharroo, Delft Institute for
# Earth-Oriented Space Research, Delft University of Technology, Delft,
# Netherlands. E-mail: remko.scharroo@lr.tudelft.nl
#
########################################################################
#FC	= f77
#CC	= cc
IFLAGS	= -I.
FFLAGS	= -O $(IFLAGS)
CFLAGS	= -O $(IFLAGS)
CPP	= /lib/cpp -P

########################################################################
SUBS	= lgrid.o gridrd.o i2swap.o fastio.o

.SUFFIXES:	.F

all:	lgrid

########################################################################
$(SUBS):	sysdep.inc sysdep.h

lgrid:	$(SUBS)
	$(FC) -o $@ $(SUBS)

.F.o:
	$(CPP) $(IFLAGS) $< > F$*.f
	$(FC) $(FFLAGS) -c F$*.f
	mv F$*.o $@
	rm -f F$*.f

sysdep.h sysdep.inc:        sysdep
	./sysdep

sysdep:       underscore.c sysdep.f
	$(CC) -c underscore.c
	$(FC) $(FFLAGS) sysdep.f -o $@ underscore.o
	rm -f underscore.o sysdep.o

test:	lgrid test.grd test.sed test.asc
	lgrid test.grd | sed -f test.sed > test.tmp
	diff test.tmp test.asc
	@echo "****************************************"
	@echo "*** Files are identical. Test passed ***"
	@echo "****************************************"

clean:
	-rm -f *.o sysdep sysdep.inc sysdep.h sysdep.mk test.tmp

spotless:	clean
	-rm -f lgrid

tar:	spotless
	( cd .. ; tar -cvhf - lgrid | gzip > lgrid.tar.gz )
