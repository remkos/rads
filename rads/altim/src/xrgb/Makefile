########################################################################
#	Makefile for XRGB
#
# 15-Feb-1998 - Multiplatform version
#
########################################################################

include $(ALTIM)/include/config.mk

########################################################################

SUBS	= addhls.o addrgb.o boxes.o \
	cmchange.o cmcopy.o cmlinear.o cmload.o cmsave.o \
	colupdt.o menu.o tag.o untag.o whatindx.o xrgb.o
EXEC	= $(BIN)/xrgb

########################################################################

all: $(EXEC)

clean:
	$(RM) -f $(EXEC) *.o

$(SUBS)	: xrgb.inc

$(BIN)/xrgb:	$(SUBS)
	$(FC) -o $@ $(SUBS) $(PMPLOT)
	strip $@
