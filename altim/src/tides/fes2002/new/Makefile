#######################################################################
# Makefile for the FES prediction software.
# Version 6.1 
# Major changes: June, 3rd 2003
# Full revision by Fabien LEFEVRE and Frédéric BRIOL
#
# This software have been tested on SunOS 5.6
# and Red Hat Linux release 7.1 (Seawolf)
# It is provided without any guarantees.
#
# For bug reports, please contact :
# ---------------------------------
# Fabien LEFEVRE 
# 
# CLS
# http://www.cls.fr
# Direction Océanographie Spatiale
# 8-10, rue Hermès - Parc Technologique du Canal
# 31526 Ramonville Saint-Agne cedex - France
# Tel: +33 (0)5 61 39 37 45 Fax: +33 (0)5 61 39 37 82
# e-mail: Fabien.Lefevre@cls.fr
# 
#######################################################################

#-- If operating system is Linux
CC		= gcc
FC		= g77
FFLAGS		= -O
CFLAGS		= -O

#-- If operating system is Sun Solaris
# CC 		= cc
# FC		= f77
# FFLAGS	= -O
# CFLAGS	= -O

#-- If operating system is SGI/IRIX
# CC 		= cc
# FC		= f77
# FFLAGS	= -O
# CFLAGS	= -O

#-- If operating system is HP/UX
# CC 		= cc
# FC		= f77
# FFLAGS	= -O
# CFLAGS	= -O -DNOUNDERSCOREAFTER

#-- If operating system is IBM/AIX
# CC 		= cc
# FC		= f77
# FFLAGS	= -O
# CFLAGS	= -O -DNOUNDERSCOREAFTER

#######################################################################
# Nothing to change in the following...

#----------------------------------------------------------------------------------
.SUFFIXES: .o .f .c

.f.o:
	$(FC) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<
#----------------------------------------------------------------------------------

BINARIES =	testfes	fcc_testfes

OBJ_C = 	fescore.o festide.o testfes.o
OBJ_F =		fescore.o festide.o fcc_testfes.o

SYS_LIB =	-lm

#-- all
all : $(BINARIES)

#-- clean
clean :
	$(RM) $(BINARIES) $(OBJ_C) $(OBJ_F)

#-- dependencies
testfes.c: 	fescore.h festide.h
fcc_testfes.f: 	fescore.h festide.h

testfes :	$(OBJ_C)
	$(CC) $(CFLAGS) $(OBJ_C) $(SYS_LIB) -o $@

fcc_testfes :	$(OBJ_F)
	$(FC) $(FFLAGS) $(OBJ_F)  $(SYS_LIB) -o $@ 
#######################################################################
