#
# Makefile for IRI95
#
include $(ALTIM)/include/config.mk
FFLAGS	+=-I$(NETCDFHOME)/include

SUBS	= iri95tec.o irifun.o irisub.o ../../rssubs/checkenv.o ../../rssubs/freeunit.o ../../rssubs/mdate.o
FFLAGS	:= $(subst -fimplicit-none,,$(FFLAGS))

all:	$(IRI95)

clean:
	$(RM) -f $(IRI95) *.o

$(IRI95):	$(SUBS)
	rm -f $@
	$(AR) $(ARFLAGS) $@ $(SUBS)
	$(RANLIB) $@
