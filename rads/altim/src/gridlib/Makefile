#####################################################################
#
# Makefile to generate the Grid library.
#
# 14-Feb-1998 - Multiplatform version
#
#####################################################################

include $(ALTIM)/include/config.mk

#####################################################################

all:	$(GRID)

clean:
		$(RM) -f $(GRID) *.o

#####################################################################

FSUBS1	= dyntopo.f griddx.f griddy.f grillu.f
FSUBS2	= ../rssubs/i2swap.f ../rssubs/i4swap.f ../rssubs/isnan.f \
	../rssubs/fin.f ../rssubs/lnblnk.f ../rssubs/ltlend.f \
	conic.f grest1.f grest2.f grintp.f grlcom.f grltra.f grsexp.f grsqrt.f
FSUBS3	= gr1int.F gr2int.F gr3int.F grwint.F gridrd.F gridwr.F sgridrd.F \
	grstat.F gridcog.F
FSUBS4	= gridlib.f90 masklib.f90
CSUBS1	= bitmask.c
CSUBS2	= ../rssubs/fastio.c ../rssubs/mallocf.c ../rssubs/memloc.c
SUBS	= $(FSUBS1:.f=.o) $(FSUBS2:.f=.o) $(FSUBS4:.f90=.o) \
	$(FSUBS3:.F=4.o) $(FSUBS3:.F=8.o) $(CSUBS1:.c=.o) $(CSUBS2:.c=.o)
TARFILES	= *.f *.F Makefile*
MAN	= $(DOCDIR)/grid.pdf
TEXSUBS	= $(FSUBS1:.f=.tex) $(FSUBS3:.F=4.tex) $(FSUBS3:.F=8.tex) $(FSUBS4:.f90=.tex)

# How to make the gridwr and gridrd subroutines from the .F versions
#
.SUFFIXES:	.f90 .f .c .F 4.o 8.o .o 4.tex 8.tex
.F4.o:
		$(CPP) $< > F$(<:.F=.f)
		$(FC) $(FFLAGS) -c F$(<:.F=.f) -o $@
		rm -f F$(<:.F=.f)
.F8.o:
		$(CPP) -DDOUBLE $< > F$(<:.F=.f)
		$(FC) $(FFLAGS) -c F$(<:.F=.f) -o $@
		rm -f F$(<:.F=.f)

.F4.tex:
		$(CPP) $< > F$(<:.F=.f)
		maketex F$(<:.F=.f) > $@
		rm -f F$(<:.F=.f)

.F8.tex:
		$(CPP) -DDOUBLE $< > F$(<:.F=.f)
		maketex F$(<:.F=.f) > $@
		rm -f F$(<:.F=.f)

# Dependancies
#
gridbint.o gridbinq.o gridbuff_deos.o gridbuff_nc gridbspl.o \
maskbint.o maskbuff.o:	$(ALTIM)/include/nan.inc gridbuff.inc

# Specific make rules
#
gridbuff_nc.o:	gridbuff_nc.f
		$(FC) $(FFLAGS) -I$(NETCDFHOME)/include gridbuff_nc.f -c -o $@

gridlib.o:	gridlib.f90
		$(FC) $(FFLAGS) -I$(NETCDFHOME)/include gridlib.f90 -c
		mv -f gridlib.mod $(ALTIM)/include

masklib.o:	masklib.f90
		$(FC) $(FFLAGS) -I$(NETCDFHOME)/include masklib.f90 -c
		mv -f masklib.mod $(ALTIM)/include

bitmask.o:      bitmask.h machine.h bitmask.c
		$(CC) -c -O $(CFLAGS) $(IFLAGS) $(CPLUS) bitmask.c

# How to make manual
#
man:	$(MAN)

$(MAN):	grid.tex $(TEXSUBS)
	pdflatex grid
	mv grid.pdf $@
	rm -f grid.{dvi,log}

# The rest of this Makefile describes how to make the grid library.
#
$(GRID):	$(SUBS)
		rm -f $@
		$(AR) $(ARFLAGS) $@ $(SUBS)
		$(RANLIB) $@

tar:		gridlib.tar.gz
gridlib.tar.gz:	$(TARFILES)
		tar -cvf - $(TARFILES) | gzip > $@

# How to make everything clean
#
clean:
		-rm -f *.o *.bak *~ *.swp

cleanup:	clean
		-rm -f $(GRID) gridlib.tar.gz
