########################################################################
#
#             Makefile for NUMLIB
#
# 15-Feb-1998 - Multiplatform version
#
########################################################################
include $(ALTIM)/include/config.mk

NUMLIB	= $(ALTIM)/lib/numlib.a
FFLAGS	= -O

########################################################################

#FSUBS	= bessel.f chepol.f cheree.f chetra.f chever.f chevet.f chevre.f diamst.f diatsm.f eigcon.f eighes.f eigjac.f eigjad.f errorf.f funafg.f gammaf.f herpol.f herree.f hertra.f intexp.f intga1.f intgau.f invcho.f inverf.f invgau.f invsyd.f kartri.f karweb.f klefun.f klepol.f lagpol.f lagree.f lagtra.f linban.f linch1.f lincho.f lingau.f linmgs.f linsy1.f linsyd.f lintri.f matban.f matbas.f matbmv.f matcho.f matdet.f matmgs.f matmmm.f matmmv.f matmst.f matsmm.f matsmv.f matssm.f matsy1.f matsyd.f matsyr.f mattmv.f mattsm.f mingrf.f minmul.f minrrf.f nuldri.f nultwe.f polbar.f poldel.f polhor.f polsub.f polver.f runeen.f runner.f runnor.f spicub.f spicun.f zerfun.f zerpol.f
FSUBS	= matsy1.f matcho.f mattmv.f linch1.f matmmv.f invcho.f matsmv.f

all: $(NUMLIB)

clean:
	$(RM) -f $(NUMLIB) *.o

########################################################################
#	Library
########################################################################

$(NUMLIB):	$(FSUBS)
		make $(?:.f=.o)
		$(AR) $(ARFLAGS) $@ $(?:.f=.o)
		$(RANLIB) $@
		rm -f $(?:.f=.o)
