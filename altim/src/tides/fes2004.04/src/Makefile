# ######################################################################
#  Makefile FES prediction software.
# 
#  File      : Makefile
#  Developer : CLS
#  Version   : 1.1
#  Date      : 6 oct 2004
#  
#  This software have been tested on Linux platform
#  It is provided without any guarantees.
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
#        developed by Jean-Marc MOLINES and Florent LYARD
# #####################################################################

# --------------------------------------------------------------------
# TO BE MODIFIED IF NEEDED !!!
# 
# Select one and only one of the three following options of compilation:
# 1) load ascii files (long initialisation, few hard drive space)
#
# 2) load binary files (need more hard drive space, but faster execution)
# 3) load binary files + direct access on the binary files
#     (need more hard drive space, but fasvery fast execution!!!) 
# Comment the other option with a # at the beginning of the line
#
# 1) load ascii files
CFLAGS=-g
#
# 2) load binary files
# CFLAGS=-O -DBINARY
#
# 3) load binary files + direct accesss on the binary file
# CFLAGS=-O -DBINARY -DIO
#
#
#-- If operating system is Linux
CC		= gcc
FC		= pgf77
#
#-- If operating system is Sun Solaris
# CC 		= cc
# FC		= f77
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
# ascii
#
ascii : lib test check

#
# all
#
binary : lib bin test check

#
# clean
#
clean :
	$(RM) $(LIB_FES) $(LIB_OBJS) *.o $(BINARIES) output_file.txt
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
testfes : testfes.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS) -lfes -lm	

fcc_testfes :  fcc_fes.o fcc_testfes.o
	$(FC) -o $@ $(CFLAGS) $^ $(LDFLAGS) -lfes -lm	

#
# ascii2bin
#
ascii2bin.o: ascii2bin.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $<

ascii2bin: ascii2bin.o
	$(CC) -o $@ $(CFLAGS) $^ -lm	

DATA	= ../data

#
# converts ASCII to BIN
#
bin:	ascii2bin		\
	$(DATA)/2N2_drfes2004.bin	\
	$(DATA)/2N2_fes2004.bin	\
	$(DATA)/K1_drfes2004.bin	\
	$(DATA)/K1_fes2004.bin	\
	$(DATA)/K2_drfes2004.bin	\
	$(DATA)/K2_fes2004.bin	\
	$(DATA)/M2_drfes2004.bin	\
	$(DATA)/M2_fes2004.bin	\
	$(DATA)/N2_drfes2004.bin	\
	$(DATA)/N2_fes2004.bin	\
	$(DATA)/O1_drfes2004.bin	\
	$(DATA)/O1_fes2004.bin	\
	$(DATA)/P1_drfes2004.bin	\
	$(DATA)/P1_fes2004.bin	\
	$(DATA)/Q1_drfes2004.bin	\
	$(DATA)/Q1_fes2004.bin	\
	$(DATA)/S2_drfes2004.bin	\
	$(DATA)/S2_fes2004.bin


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
