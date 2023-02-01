#
# Makefile for IRI2007
#
include $(ALTIM)/include/config.mk

SUBS	= iri2007tec.o irifun.o irisub.o igrf.o ../../rssubs/checkenv.o ../../rssubs/freeunit.o ../../rssubs/mdate.o ../../rssubs/isnan.o
FFLAGS	:= $(subst -fimplicit-none,,$(FFLAGS))

all:	$(IRI2007)

clean:
	$(RM) -f $(IRI2007) *.o

$(IRI2007):	$(SUBS)
	rm -f $@
	$(AR) $(ARFLAGS) $@ $(SUBS)
	$(RANLIB) $@
