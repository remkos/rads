########################################################################
#
#	Makefile for RSSUBS
#
# 15-Feb-1998 - Multiplatform version
#
########################################################################
ALTIM=/Users/davidtrossman/Documents/Code/rads/altim
HOSTTYPE=x86_64
#HOSTTYPE=i386_gfortran
#NETCDFHOME=/opt/local
NETCDFHOME=/opt/local/lib

include $(ALTIM)/include/config.mk
CPPFLAGS=

########################################################################
# NO ADJUSTMENTS should be needed beyond this point.
########################################################################

########################################################################
#	Definition of object generation
########################################################################
.SUFFIXES:	.F .tex

########################################################################
#	Source definitions
########################################################################
FSUBS	= altbias.f anmtot.f anttom.f bubble.f \
	carpol.f checkenv.f chartrans.f chrdat.f chrloc.f covar.f \
	datearg.f dhellips.f dqsort.f dyntopo.f \
	earth.f elevec.f epharg.f ephpos.f f1f2.f fin.f freeunit.f \
	gauss1d.f geocen.f geodet.f geoutm.f geoxyz.f \
	getorb.f globpres.f grdate.f helmert1.f helmert2.f \
	i2swap.f i4swap.f i8swap.f intab2.f intab8.f inter8.f iqsort.f isnan.f \
	j2000.f land.f listargs.f lnblnk.f loadcoord.f lowercase.f ltlend.f \
	matsy1.f mcw.f mdate.f mjdate.f nff.f noise.f \
	nuvel1a.f odrinfo.f polcar.f quaint.f regres.f rotate.f round.f rqsort.f \
	scaprd.f sec85.f sfdist.f sland.f statbar.f statinfo.f statis.f \
	strf1985f.f thetop.f uppercase.f vecele.f vecprd.f vecnrm.f vnorm.f \
	xyzgeo.f ymdhms.f
CSUBS	= fastio.c fdate.c mallocf.c memloc.c strfepoch.c strpepoch.c
F90SUBS = nf_subs.f90 solar_subs.f90 tpj_subs.f90 meteo_subs.f90
OBJS	= $(FSUBS:.f=.o) $(CSUBS:.c=.o) $(F90SUBS:.f90=.o)
TEXSUBS	= $(OBJS:.o=.tex)
TARFILES = Makefile* *.[fF] *.inc *.c rssubs.tex
MAN	= $(DOCDIR)/rssubs.pdf
EXEC	= $(BIN)/mdate $(BIN)/xyzgeo $(BIN)/geoxyz $(BIN)/convdate \
	$(BIN)/altbias $(BIN)/getorb

########################################################################
#	What to make?
########################################################################

all:	lib exec
lib:	$(RSSUBS)
exec:	$(EXEC)
clean:
		$(RM) -f *.o $(RSSUBS) $(EXEC)

########################################################################
#	Include files
########################################################################

anmtot.o vecele.o: $(ALTIM)/include/math.inc
covar.o: covar.inc
geoutm.o covar.inc: $(ALTIM)/include/GRS80.inc
earth.o geocen.o geodet.o xyzgeo.o: earth.inc

########################################################################
#	Object generation
########################################################################

iqsort.f:	qsort.F
	$(CPP) qsort.F > iqsort.f

dqsort.f:	qsort.F
	$(CPP) -DDOUBLE qsort.F > dqsort.f

rqsort.f:	qsort.F
	$(CPP) -DSINGLE qsort.F > rqsort.f

nff.o:		nff.f
	$(FC) $(FFLAGS) -I$(NETCDFHOME)/include nff.f -c

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NETCDFHOME)/include $*.f90 -c
	if test -f $*.mod; then mv -f $*.mod $(ALTIM)/include; fi

$(RSSUBS):	$(OBJS)
	rm -f $@
	mkdir -p $(LIBDIR)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@

$(BIN)/mdate:	mdate_main.o $(RSSUBS)
	$(FC) -o $@ mdate_main.o $(RSSUBS)
	strip $@

$(BIN)/xyzgeo:	xyzgeo_main.o $(RSSUBS)
	$(FC) -o $@ xyzgeo_main.o $(RSSUBS)
	strip $@

$(BIN)/geoxyz:	geoxyz_main.o $(RSSUBS)
	$(FC) -o $@ geoxyz_main.o $(RSSUBS)
	strip $@

$(BIN)/convdate:	convdate.o $(RSSUBS)
	$(FC) -o $@ convdate.o $(RSSUBS)
	strip $@

$(BIN)/altbias:	altbias_main.o $(RSSUBS)
	$(FC) -o $@ altbias_main.o $(RSSUBS)
	strip $@

$(BIN)/getorb:	getorb_main.o $(RSSUBS)
	$(FC) -o $@ getorb_main.o $(RSSUBS)
	strip $@

helmert: helmert.o
	$(FC) -o $@ helmert.o $(RSSUBS) $(MATHLIB)

tar:	rssubs.tar.gz
rssubs.tar.gz:	$(TARFILES)
	tar -cvf - $(TARFILES) | gzip > $@

man:	$(MAN)

$(MAN):	rssubs.tex $(TEXSUBS)
	pdflatex rssubs
	mv rssubs.pdf $@
	rm -f rssubs.{dvi,log}

clean:
	-rm -f *.o $(TEXSUBS) syscheck sysinfo *.bak *~ core

cleanup:	clean
	rm -f $(RSSUBS) $(EXEC) rssubs.tar.gz
