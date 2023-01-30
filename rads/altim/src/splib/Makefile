#####################################################################
#
# Makefile to create SPLIB
#
# 14-Feb-1998 - Multiplatform version
#
#####################################################################
# System definitions

include $(ALTIM)/include/config.mk

#####################################################################
# Defining the subroutines
#
FSUBS	= spcdem.f spcdemn.f spcomp.f spdftc.f spdftr.f spectr.f \
	spfamp.f spfftc.f spfftr.f spfper.f sprand.f statis.f

#####################################################################
# What should be zipped
#
TARFILES= *.f *.tex *.aux *.toc *.log *.dvi Makefile

#####################################################################
# What are the tex files
#
TEXSUBS	= $(FSUBS:.f=.tex)
MAN	= $(DOCDIR)/splib.pdf

#####################################################################
# What to do on `make all'
#
all:	$(SPLIB)

clean:
	$(RM) -f $(SPLIB) *.o

#####################################################################
# Rest of the defines
#
.SUFFIXES:	.tex

$(SPLIB):	$(FSUBS)
		$(FC) $(FFLAGS) -c $?
		$(AR) $(ARFLAGS) $@ $(?:.f=.o)
		$(RANLIB) $@
		rm -f $(?:.f=.o)

man:		$(MAN)

$(MAN):		splib.tex $(TEXSUBS)
		pdflatex splib
		mv splib.pdf $@
		rm -f splib.{dvi,log}

tar:		splib.tar.gz
splib.tar.gz:	$(TARFILES)
		tar -cvf $@ $(TARFILES)
