########################################################################
#
#	Makefile for SYSDEP
#
# Creates system dependent files: sysdep.h, sysdep.inc, sysdep.mk
# and support programs
#
########################################################################

ALTIM	= /user/altim
BIN	= $(ALTIM)/bin
FFLAGS	= -O
CFLAGS	= -O
SYSDEPS	= sysdep.h sysdep.mk
OBJS	= sysdep.o underscore.o

all:		sysdep.h

$(SYSDEPS):	sysdep
		./sysdep

sysdep:		$(OBJS)
		$(FC) -o $@ $(OBJS)
		strip $@

install:
		cp -p $(SYSDEPS) $(ALTIM)/include

clean:
		rm -f $(SYSDEPS) *.o sysdep
