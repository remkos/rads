########################################################################
#
#	Makefile for TOOLS
#
# Create support programs
#
########################################################################

include $(ALTIM)/include/config.mk

########################################################################
# NO ADJUSTMENTS should be needed beyond this point.
########################################################################

EXEC	= $(BIN)/maketex $(BIN)/makehtml $(BIN)/pgpline

all:	$(EXEC)

clean:
		$(RM) -f $(EXEC) *.o

$(BIN)/maketex:	maketex.o lnblnk.o
		$(FC) -o $@ maketex.o lnblnk.o
		strip $@

$(BIN)/makehtml:	$(BIN)/maketex
		rm -f $@
		ln -s maketex $@

$(BIN)/pgpline:	pgpline.o lnblnk.o
		$(FC) -o $@ pgpline.o lnblnk.o
		strip $@
