#######################################################################
# Makefile for the FES prediction software.
# Version 3.1
# Major changes: full revision by Fabien LEFEVRE
#
# Files are distributed as ASCII files, but can be converted to binary
# files by the program ascii2bin. Doing so reduces both the required disk
# space and system time for reading the files.
# Depending on your choice, modify the CPPFLAfestide.FG in the Makefile
# -DBIN_FILE : to compile to use binary file
# -DADM_COMP : to compile to use computed damittance wave (not gridded wave)
#              => the soft will run faster as only 8 grids will be read
#
# This software have been tested on SUN Solaris 2.7
# It is provided without any guarantees.
#
#   For bug reports, please contact :
#   ---------------------------------
#   Fabien LEFEVRE 
# 
#   CLS
#   http://www.cls.fr
#   Direction Océanographie Spatiale
#   8-10, rue Hermès - Parc Technologique du Canal
#   31526 Ramonville Saint-Agne cedex - France
#   Tel: +33 (0)5 61 39 37 45 Fax: +33 (0)5 61 39 37 82
#   e-mail: Fabien.Lefevre@cls.fr
# 
#   NOTE: This software is based on the former versions
#         developed by Jean-Marc MOLINES and Florent LYARD
#######################################################################

#-- Optimized flags
FFLAGS= -O $(CPPFLAG) 
CFLAGS= -O $(CPPFLAG) 

#-- Debug flags
#FFLAGS= -g -u -C $(CPPFLAG) 
#CFLAGS= -g $(CPPFLAG) 

FC=f77
CC=cc

EXECDIR=../exe

#######################################################################

.SUFFIXES: .o .f .c

.f.o:
	$(FC) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(FFLAGS) $<

#######################################################################

all: fes testfes testfes_c ascii2bin libfestide

fes: festide.o lpeqmt.o fes.o
	$(FC) $(FFLAGS) fes.o festide.o lpeqmt.o -o $(EXECDIR)/fes

testfes: festide.o lpeqmt.o testfes.o
	$(FC) $(FFLAGS) testfes.o festide.o lpeqmt.o -o $(EXECDIR)/testfes

testfes_c: festide.o lpeqmt.o testfes_c.o
	$(FC) $(FFLAGS) testfes_c.o festide.o lpeqmt.o -o $(EXECDIR)/testfes_c

ascii2bin: ascii2bin.o festide.o 
	$(FC) $(FFLAGS) ascii2bin.o festide.o lpeqmt.o -o $(EXECDIR)/ascii2bin

libfestide: festide.o lpeqmt.o 
	ar -r libfestide.a festide.o lpeqmt.o 

clean:
	\rm -f $(OBJS) *.o core *.% *.a *.*~
	
# dependencies
festide.o : common.h
fes.o     : festide.o
testfes.o : festide.o
