########################################################################
#
#       Makefile for INGRID
#
# 15-Feb-1998 - Multiplatform version
#
########################################################################

include $(ALTIM)/include/config.mk

########################################################################
# NO ADJUSTMENTS should be needed beyond this point.
########################################################################

SUBS1	= cparea.o ingrid.o ireadf.o iwritef.o menuh.o newbuf.o noclob.o \
	oldbuf.o poly.o repack.o twobuf.o
SUBS2	= getrange.o multin.o noneg.o rot90.o subarr.o wieber.o
SUBS	= $(SUBS1) $(SUBS2)
EXEC	= $(BIN)/ingrid

########################################################################

all:	$(EXEC)
clean:
		$(RM) -f $(EXEC) *.o

$(SUBS1):	ingrid.inc

#ingrid.o:	ingrid.f
#	$(FC) $(FFLAGS) -NQ20000 -c ingrid.f

$(BIN)/ingrid:	$(SUBS) $(GRID) $(RSSUBS)
	$(FC) -o $@ $(SUBS) $(GRID) $(RSSUBS)
	strip $@
