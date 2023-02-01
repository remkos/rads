#
# Makefile for ionospheric models
#
ALTIM=/Users/davidtrossman/Documents/Code/rads/altim
HOSTTYPE=x86_64
#HOSTTYPE=i386_gfortran
#NETCDFHOME=/opt/local
NETCDFHOME=/opt/local/lib
include $(ALTIM)/include/config.mk

FFLAGS	+=-I$(NETCDFHOME)/include

SUBS	= iritec.o gimread.o gimtec.o gimsubs.o nic07tec.o nicsubs.o \
	../rssubs/lnblnk.o ../rssubs/mdate.o \
	../rssubs/freeunit.o ../rssubs/checkenv.o
MODULES	= ../../include/nicsubs.mod ../../include/gimsubs.mod

%.o:	%.mod

all:	$(IONO) $(MODULES)

clean:
		$(RM) -f $(IONO) $(MODULES) *.o *.mod

$(IONO):	$(SUBS)
	rm -f $@
	$(AR) $(ARFLAGS) $@ $(SUBS)
	$(RANLIB) $@

gimtest:	gimtest.o $(IONO)
	$(FC) -o $@ gimtest.o $(IONO) $(NETCDF)

$(MODULES):	nicsubs.o gimsubs.o
	cp nicsubs.mod gimsubs.mod ../../include
