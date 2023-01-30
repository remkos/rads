########################################################################
#
#       Makefile for GEODYN Tools
#
# 18-Feb-1998 - Multiplatform version
#
########################################################################

include $(ALTIM)/include/config.mk
FFLAGS=-O

########################################################################
# NO ADJUSTMENTS should be needed beyond this point.
########################################################################

all:		$(BIN)/tle2elem

clean:
	$(RM) -f $(BIN)/tle2elem *.o

########################################################################

$(BIN)/tle2elem: tle2elem.o tymd.o uneof.o ymdt.o $(RSSUBS)
		$(FC) -o $@ tle2elem.o tymd.o uneof.o ymdt.o $(RSSUBS)
		strip $@ 

