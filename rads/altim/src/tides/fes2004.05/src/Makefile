# ######################################################################
#  Makefile FES prediction software.
# 
#  File      : Makefile
#  Developer : CLS
#  Version   : 1.3
#  Date      : 19 January 2005
#  
#  This software have been tested on Linux platform
#  It is provided without any warranty.
# 
#  For bug reports, please contact :
#  ---------------------------------
#  Fabien LEFEVRE 
# 
#  CLS
#  http://www.cls.fr
#  Direction Océanographie Spatiale
#  8-10, rue Hermès - Parc Technologique du Canal
#  31526 Ramonville Saint-Agne cedex - France
#  Tel: +33 (0)5 61 39 37 45 Fax: +33 (0)5 61 39 37 82
#  e-mail: Fabien.Lefevre@cls.fr
# 
#  NOTE: This software is based on the former versions
#        developed by CNRS (Jean-Marc MOLINES and
#        Florent LYARD among others)
# #####################################################################

# --------------------------------------------------------------------
# TO BE MODIFIED IF NEEDED !!!
# 
CFLAGS=-O
#CFLAGS=-g
#
#-- If operating system is Linux
CC		= gcc
FC		= pgf77
#
#-- If operating system is Sun Solaris
#CC 		= cc
#FC		= f77
#
# End of TO BE MODIFIED IF NEEDED !!!
# --------------------------------------------------------------------

INCLUDE = -I../include
LIB_DIR = ../lib
LIB_FES = $(LIB_DIR)/libfes.a
LDFLAGS = -L$(LIB_DIR)


.SUFFIXES: .o .f .c

.c.a:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $*.o
	$(AR) r $@ $*.o
	$(RM) $*.o

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $<

.f.o:
	$(FC) -c $(CFLAGS) $(INCLUDE) $<


LIB_OBJS =	interp.o	\
		error.o		\
		grid.o		\
		prediction.o	\
		fes.o
		
BINARIES = ascii2bin testfes fcc_testfes

#
# all
#
all : bin lib bin test check

#
# ascii
#
ascii : lib test check

#
# clean
#
clean :
	$(RM) $(LIB_FES) $(LIB_OBJS) *.o $(BINARIES) fes.output.fes2004.test
#
# lib
#
lib : $(LIB_FES)

#
# test
#
test : lib $(BINARIES)

#
# check
#
check : test
	./testfes > fes.output.fes2004.test
	cmp fes.output.fes2004 fes.output.fes2004.test


#
# library depends
#
$(LIB_FES)(interp.o) :		interp.c
$(LIB_FES)(error.o) :		error.c
$(LIB_FES)(grid.o) :		grid.c
$(LIB_FES)(prediction.o) :	prediction.c
$(LIB_FES)(lpeqmt.o) :		lpeqmt.c
$(LIB_FES)(fes.o) :		fes.c

$(LIB_FES) :	$(LIB_FES)(interp.o)
$(LIB_FES) :	$(LIB_FES)(error.o)
$(LIB_FES) :	$(LIB_FES)(grid.o)
$(LIB_FES) :	$(LIB_FES)(prediction.o)
$(LIB_FES) :	$(LIB_FES)(fes.o)

#
# testFes
#
testfes : testfes.o $(LIB_FES)
	$(CC) -o $@ $(CFLAGS) testfes.o $(LDFLAGS) -lfes -lm	

fcc_testfes : fcc_fes.o fcc_testfes.o $(LIB_FES)
	$(FC) -o $@ $(CFLAGS) fcc_fes.o fcc_testfes.o $(LDFLAGS) -lfes -lm	

#
# ascii2bin
#
ascii2bin : ascii2bin.o
	$(CC) -o $@ $(CFLAGS) ascii2bin.o -lm 

DATA	= ../data

#
# converts ASCII to BIN
#
bin:	ascii2bin			\
	$(DATA)/2N2_drfes2004.bin	\
	$(DATA)/2N2_fes2004.bin		\
	$(DATA)/K1_drfes2004.bin	\
	$(DATA)/K1_fes2004.bin		\
	$(DATA)/K2_drfes2004.bin	\
	$(DATA)/K2_fes2004.bin		\
	$(DATA)/M2_drfes2004.bin	\
	$(DATA)/M2_fes2004.bin		\
	$(DATA)/N2_drfes2004.bin	\
	$(DATA)/N2_fes2004.bin		\
	$(DATA)/O1_drfes2004.bin	\
	$(DATA)/O1_fes2004.bin		\
	$(DATA)/P1_drfes2004.bin	\
	$(DATA)/P1_fes2004.bin		\
	$(DATA)/Q1_drfes2004.bin	\
	$(DATA)/Q1_fes2004.bin		\
	$(DATA)/S2_drfes2004.bin	\
	$(DATA)/S2_fes2004.bin		\
	$(DATA)/M4_fes2004.bin		\
	$(DATA)/S1_fes2004.bin		\
	$(DATA)/Mf_fes2004.bin		\
	$(DATA)/Mm_fes2004.bin		\
	$(DATA)/Mtm_fes2004.bin 	\
	$(DATA)/MSqm_fes2004.bin


#
# binary depends
#
$(DATA)/2N2_drfes2004.bin : $(DATA)/2N2_drfes2004.asc 
	ascii2bin $(DATA)/2N2_drfes2004.asc $(DATA)/2N2_drfes2004.bin 0
$(DATA)/2N2_fes2004.bin : $(DATA)/2N2_fes2004.asc 
	ascii2bin $(DATA)/2N2_fes2004.asc $(DATA)/2N2_fes2004.bin 0
$(DATA)/K1_drfes2004.bin : $(DATA)/K1_drfes2004.asc 
	ascii2bin $(DATA)/K1_drfes2004.asc $(DATA)/K1_drfes2004.bin 0
$(DATA)/K1_fes2004.bin : $(DATA)/K1_fes2004.asc 
	ascii2bin $(DATA)/K1_fes2004.asc $(DATA)/K1_fes2004.bin 0
$(DATA)/K2_drfes2004.bin : $(DATA)/K2_drfes2004.asc 
	ascii2bin $(DATA)/K2_drfes2004.asc $(DATA)/K2_drfes2004.bin 0
$(DATA)/K2_fes2004.bin : $(DATA)/K2_fes2004.asc 
	ascii2bin $(DATA)/K2_fes2004.asc $(DATA)/K2_fes2004.bin 0
$(DATA)/M2_drfes2004.bin : $(DATA)/M2_drfes2004.asc 
	ascii2bin $(DATA)/M2_drfes2004.asc $(DATA)/M2_drfes2004.bin 0
$(DATA)/M2_fes2004.bin : $(DATA)/M2_fes2004.asc 
	ascii2bin $(DATA)/M2_fes2004.asc $(DATA)/M2_fes2004.bin 0
$(DATA)/N2_drfes2004.bin : $(DATA)/N2_drfes2004.asc 
	ascii2bin $(DATA)/N2_drfes2004.asc $(DATA)/N2_drfes2004.bin 0
$(DATA)/N2_fes2004.bin : $(DATA)/N2_fes2004.asc 
	ascii2bin $(DATA)/N2_fes2004.asc $(DATA)/N2_fes2004.bin 0
$(DATA)/O1_drfes2004.bin : $(DATA)/O1_drfes2004.asc 
	ascii2bin $(DATA)/O1_drfes2004.asc $(DATA)/O1_drfes2004.bin 0
$(DATA)/O1_fes2004.bin : $(DATA)/O1_fes2004.asc 
	ascii2bin $(DATA)/O1_fes2004.asc $(DATA)/O1_fes2004.bin 0
$(DATA)/P1_drfes2004.bin : $(DATA)/P1_drfes2004.asc 
	ascii2bin $(DATA)/P1_drfes2004.asc $(DATA)/P1_drfes2004.bin 0
$(DATA)/P1_fes2004.bin : $(DATA)/P1_fes2004.asc 
	ascii2bin $(DATA)/P1_fes2004.asc $(DATA)/P1_fes2004.bin 0
$(DATA)/Q1_drfes2004.bin : $(DATA)/Q1_drfes2004.asc 
	ascii2bin $(DATA)/Q1_drfes2004.asc $(DATA)/Q1_drfes2004.bin 0
$(DATA)/Q1_fes2004.bin : $(DATA)/Q1_fes2004.asc 
	ascii2bin $(DATA)/Q1_fes2004.asc $(DATA)/Q1_fes2004.bin 0
$(DATA)/S2_drfes2004.bin : $(DATA)/S2_drfes2004.asc 
	ascii2bin $(DATA)/S2_drfes2004.asc $(DATA)/S2_drfes2004.bin 0
$(DATA)/S2_fes2004.bin : $(DATA)/S2_fes2004.asc 
	ascii2bin $(DATA)/S2_fes2004.asc $(DATA)/S2_fes2004.bin 0
$(DATA)/M4_fes2004.bin : $(DATA)/M4_fes2004.asc 
	ascii2bin $(DATA)/M4_fes2004.asc $(DATA)/M4_fes2004.bin 0
$(DATA)/Mf_fes2004.bin : $(DATA)/Mf_fes2004.asc 
	ascii2bin $(DATA)/Mf_fes2004.asc $(DATA)/Mf_fes2004.bin 0
$(DATA)/Mm_fes2004.bin : $(DATA)/Mm_fes2004.asc 
	ascii2bin $(DATA)/Mm_fes2004.asc $(DATA)/Mm_fes2004.bin 0
$(DATA)/Mtm_fes2004.bin : $(DATA)/Mtm_fes2004.asc 
	ascii2bin $(DATA)/Mtm_fes2004.asc $(DATA)/Mtm_fes2004.bin 0
$(DATA)/MSqm_fes2004.bin : $(DATA)/MSqm_fes2004.asc 
	ascii2bin $(DATA)/MSqm_fes2004.asc $(DATA)/MSqm_fes2004.bin 0
$(DATA)/S1_fes2004.bin : $(DATA)/S1_fes2004.asc 
	ascii2bin $(DATA)/S1_fes2004.asc $(DATA)/S1_fes2004.bin 0
