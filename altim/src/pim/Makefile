########################################################################
#
# Makefile for PIM
#
# 14-Feb-1998 - Multiplatform version
#
########################################################################

include $(ALTIM)/include/config.mk

########################################################################

CMDSUBS	= cmdbackg.o cmdbad.o cmdbathy.o cmdbox.o cmdcircl.o \
	cmdcoast.o cmdcog.o cmdconto.o cmdcurso.o \
	cmddata.o cmddevic.o cmddithe.o cmddynto.o cmdtrack.o \
	cmdgpstr.o cmdwaypo.o \
	cmdillum.o cmdland.o cmdmesh.o \
	cmdlegen.o cmdline.o cmdlogo.o cmdannot.o cmdprofi.o cmdproje.o \
	cmdrange.o cmdrect.o cmdsurfa.o \
	cmdsymbo.o cmdtext.o cmdtitle.o cmdvalue.o \
	cmdveloc.o cmdviewp.o cmdxgf.o cmdpie.o
SUBS	= $(CMDSUBS) \
	copy.o dumppix.o \
	gridillu.o gridinfo.o gridread.o gridstck.o \
	illumpix.o annotate.o pgclh.o pgclv.o pgtrend.o \
	pgpxcont.o pim.o pimcol.o pimcont.o pixinter.o grdinter.o pimopt.o \
	chropt.o popcmd.o prepro1.o prepro2.o prepro3.o \
	readpgm.o readpgm1.o roundoff.o setarea.o setcol.o \
	showtrnd.o strip.o surf2pix.o xgfstck.o velocity.o \
	../xover/xtflimits.o
MAN	= $(DOCDIR)/pim.pdf

########################################################################

all:	$(BIN)/pim

clean:
	$(RM) -f $(BIN)/pim *.o

########################################################################

annotate.o $(CMDSUBS) dumppix.o gridillu.o gridread.o illumpix.o \
	pgclh.o pgclv.o pim.o pimcol.o pimcont.o setarea.o setcol.o \
	pgtrend.o showtrnd.o surf2pix.o:	pim.inc
cmdprofi.o:	cmdprofi.inc
cmdconto.o pimcont.o:	pimcont.inc

pim $(BIN)/pim:	$(SUBS) $(RSSUBS) $(GRID)
	$(FC) -o $@ $(SUBS) $(PMPLOT) $(RSSUBS) $(GRID)
	strip $@

man:	$(MAN)

.f.tex:;

$(MAN):	pim.tex
	pdflatex pim
	mv pim.pdf $@
	rm -f pim.{dvi,log}

clean:
	-rm -f $(SUBS)
