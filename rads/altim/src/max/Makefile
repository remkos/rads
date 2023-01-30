########################################################################
#
#       Makefile for MAX
#
########################################################################

include $(ALTIM)/include/config.mk

########################################################################

SUBS1	= blload.o blread.o \
	daread.o decl.o inclin.o insect.o \
	max.o nameid.o quasob.o rlam.o sortrk.o \
	wrixaf.o wrixtf.o wrixxf.o xofind.o xovdet.o
SUBS2	= arglat.o blnew.o bltween.o contlon.o grouplon.o nag.o \
	toppar.o
EXEC	= $(BIN)/max

########################################################################

all:	$(EXEC)

clean:
		$(RM) -f $(EXEC) *.o

$(SUBS1): maxcom.inc

nag.o:	nag.f
	$(FC) -O -c nag.f

$(BIN)/max max:	$(RSSUBS) $(SUBS1) $(SUBS2)
	$(FC) -o $@ $(SUBS1) $(SUBS2) $(RSSUBS)
	strip $@

man:	$(DOCDIR)/xformats.pdf

$(DOCDIR)/xformats.pdf:	xformats.tex
	pdflatex xformats
	mv xformats.pdf $@
	rm -rf xformats.{dvi,log}
