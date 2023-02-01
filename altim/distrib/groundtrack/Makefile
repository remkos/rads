include $(ALTIM)/include/defaults.mk

all:		groundtrack elemplot orbelem

groundtrack:	groundtrack.o $(RSSUBS)
		$(FC) -o $@ groundtrack.o $(RSSUBS)

orbelem:	orbelem.o $(RSSUBS) $(NUMLIB)
		$(FC) -o $@ orbelem.o $(RSSUBS) $(NUMLIB)

elemplot:	elemplot.o
		$(FC) -o $@ elemplot.o $(PMPLOT)
