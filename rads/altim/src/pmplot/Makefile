#####################################################################
#
# Makefile to generate the PMPLOT library.
#
# 15-Feb-1998 - Multiplatform version
#
#####################################################################
# NOTE that the environment variables ALTIM and PGPLOT_DIR should
# have been set!

include $(ALTIM)/include/config.mk

#####################################################################
# Fill in the right Fortran compiler flags.
# -O : Optimisation is on.
# -u : Assume "implicit none"
# -I$(PGPLOT_DIR) : Get additional includes from PGPLOT_DIR
# 
FFLAGS	+= -I$(PGPLOT_DIR)

#####################################################################
# Fill in the right C compiler flags.
# -O       : Optimisation is on.
# -DUNDERSCOREAFTER : Required if Fortran subroutines are called by C with
#            underscore after subroutine name
# 
#CFLAGS	= -O

#####################################################################
# Which routines are there ?
#
OBJS1	= grclip.o grlin0.o grmker.o grsy00.o grsyxd.o
OBJS2	= pgclip.o pgnorm.o pgnumb.o pgpt.o pgptx.o pgsvpx.o
OBJS3	= pmbar.o pmbox.o pmcinv.o pmconv.o pmcpol.o pmcpoly.o pmcvec.o pmdef.o \
	pmqdef.o pmqinf.o pmqwin.o pmrnd.o pmswin.o pmwdb.o pmwindow.o \
	pmx.o pmy.o
OBJS4	= cmaprd.o grtcir.o rngcir.o satobs.o grfileio.o
OBJS5	= ../rssubs/i2swap.o ../rssubs/i4swap.o ../rssubs/lnblnk.o ../rssubs/ltlend.o
OBJS	= $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5)

#####################################################################
# Which libraries are to be made ?
#
LIB	= libpmplot.a
SHARED_LIB	= $(LIBDIR)/libpmplot.$(SHARED_EXT)
MAN	= $(DOCDIR)/pmplot.pdf


#####################################################################
# Which tex files should be created from the subroutines ?
#
.SUFFIXES:	.tex .F .ps .pdf

TEXSUBS		= $(OBJS2:.o=.tex) $(OBJS3:.o=.tex) $(OBJS4:.o=.tex)
TEXMAIN		= pmplot.tex formats.tex ranks.tex pmplot.ver

#####################################################################
# How to handle .F files ?
#
.F.o:
		$(CPP) $< > F$*.f
		$(FC) $(FFLAGS) -c F$*.f
		mv F$*.o $@
		rm -f F$*.f

#####################################################################
# Which files should be tarred
#
TARFILES	= *.f *.inc *.c *.txt $(TEXMAIN) fonts Makefile

#####################################################################
# What to do on `make all'
#
all:		$(LIB) $(SHARED_LIB) grfont.dat example symbol
man:		$(MAN)

#####################################################################
# Objects using include files
#
pmbar.o pmbox.o pmcinv.o pmconv.o pmcpoly.o pmcvec.o \
	pmdef.o pmqdef.o pmqinf.o pmwindow.o pmx.o \
	pmy.o pmrnd.o cmaprd.o pmwdb.o: pmplot.inc
pgclip.o grmker.o grclip.o          : grclip.inc
grfileio.o:	grfileio.c

#####################################################################
# This part describes how to make the PMPLOT library.
#
$(LIB):		$(PGPLOT_DIR)/libpgplot.a $(OBJS)
		cp -p $(PGPLOT_DIR)/libpgplot.a $@
		$(AR) $(ARFLAGS) $@ $(OBJS)
		$(RANLIB) $@

$(SHARED_LIB):	$(LIB)
		-rm -rf $@ tmp
		mkdir -p tmp
		(cd tmp ; ar -x ../$(LIB))
		$(SHARED_LD) -o $@ tmp/*.o $(SHARED_LIB_LIBS)
		-rm -rf tmp

#####################################################################
# To make the manual
#
example:	example.o
		$(FC) example.o -o $@ $(PMPLOT)
		strip $@
example.ps:	example
		example example.ps/ps

symbol:		symbol.o
		$(FC) symbol.o -o $@ $(PMPLOT)
		strip $@
symbol.ps:	symbol
		symbol symbol.ps/ps 750

$(MAN):		$(TEXMAIN) $(TEXSUBS) example.pdf symbol.pdf
		pdflatex pmplot
		pdflatex pmplot
		mv pmplot.pdf $@
		rm -f pmplot.{dvi,log}

.ps.pdf:
		ps2raster -Tf $*.ps

#####################################################################
# To make a tar.gz-file
#
tar:		pmplot.tar.gz
pmplot.tar.gz:	$(TARFILES)
		tar -cvf- $(TARFILES) | gzip > $@


#####################################################################
# Target "grfont.dat" is the binary font file.
#
FONTS	= fonts/*.txt

pgpack: pgpack.o
	$(FC) -o $@ pgpack.o
	strip $@

grfont.dat: pgpack $(FONTS)
	rm -f grfont.dat
	cat $(FONTS) | ./pgpack

#####################################################################
# Clean up stuff
#
clean:
	-rm -f $(OBJS) $(TEXSUBS)
	-rm -f *.toc *.log *.aux

spotless:	clean
	rm -f libpmplot.a pgpack grfont.dat symbol example pmplot.dvi pmplot.ps
