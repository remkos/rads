# Makefile. GETORB. Version 2.2.1.  10 January 2006
#
# This Makefile contains the description for the generation and testing
# of 'getorb' and 'lodr', two programs to handle the Delft Orbital Data
# records (ODRs).
# - getorb: interpolate DUT orbits at requested epoch.
# - lodr  : list contents of DUT orbit files.
#
# Additionally, some tools and libraries are provided:
# - getorb.a : Library used for getorb tools
# - mdate    : A simple date conversion tool
#
# Optionally, this package provides
# - convdate : A powerfool date conversion tool
#
# Please read the instructions below carefully before installing the
# software. It has been attempted to make the software completely
# platform independent, and so we expect that you can use this Makefile
# without any major change.
#
# 1. After compiling and testing of the software, you may like to put the
#    executable in a 'bin' directory that is in your path.
#    
#BIN_DIR		= /home/altim/bin
BIN_DIR		= /usr/local/bin
#
# 2. 'make' is usually able to determine which is your Fortran compiler and
#    C compiler. However, if 'make' crashes, because it can't find the
#    compilers, specify the proper compiler names.
#    The C-preprocessor usually has to be defined.
#
#FC		= fc		# Generic
#FC		= g77		# Linux, older gcc versions
#FC             = gfortran      # Linux, gcc version 4.0 and higher
#
#CC		= cc		# Generic
#CC		= gcc		# Linux
#
# 3. Fortran and C compile flags are usually set by make. If you have
#    other preferences, change the next three lines and add the
#    prefered values. The -I. flag is required!
#
#CFLAGS		= -O
#FFLAGS		= -O -fullwarn
#LFLAGS		=
IFLAGS		= -I.
#
#    For HP Unix systems, you have to uncomment the next two lines
#
#FFLAGS		= -O +U77
#LFLAGS		= -lU77
#
#    For Sun Unix, uncomment the next line
#
#FFLAGS		= -native -fsimple -O4 -Nn5000 -I.
#
#    For SGI workstatios, uncomment the next line
#
#FFLAGS		= -O -bytereclen
#
# 4. The way to make an object library is to use "ar". Normally the default
#    is already set and it adds the contents table. Just to be sure, I add
#    here what may work on all platforms. Sometimes the table is still not
#    added and you need ranlib, but on others it is not available. When you
#    don't have ranlib, you don't need it and replace it by touch.
#
AR		= ar
ARFLAGS		= rv
RANLIB		= ranlib
#
# 5. On Windows systems output files need to by Unixised before comparison.
#    To do this right, change "touch" to "dos2unix" below.
#
DOS2UNIX	= touch
#DOS2UNIX	= dos2unix
#
# 6. This should be all that MIGHT need changing. Now type 'make' to create
#    the programs  'getorb'  and  'lodr'. If this is successful, you may like
#    to test it by running 'make test' and then commit it by typing
#    'make install'.
#    To create convdate, type 'make convdate'.
#    Afterwards, type 'make clean' to remove unnecessary remains.
#
# For any comments on the Makefile, feel free to contact the author.
#
# Remko Scharroo         remko@deos.tudelft.nl
#
# Delft Institute for Earth-Oriented Space Research
# Delft University of Technology
# Kluyverweg 1, 2629 HS, Delft, The Netherlands
###############################################################################
#
# It shouldn't be necessary to change anything in this part. Nevertheless,
# there's no harm in explaining what it is we are trying to establish here.
#
# A. Define all routines needed to create  'getorb'. Then specify the
#    names of the files that go into the tar.gz file and the executables
#    that are created.
#
OBJS		= carpol.o chrdat.o earth.o freeunit.o getorb.o \
		geocen.o geoxyz.o i2swap.o i4swap.o intab8.o inter8.o ltlend.o \
		mdate.o odrinfo.o polcar.o sec85.o fin.o datearg.o lnblnk.o

SOURCES		= getorb/index* getorb/*.[fc] getorb/*.inc \
		getorb/Makefile getorb/data getorb/docs
EXECUTABLES	= getorb lodr mdate
#
# B. Define how to make .o files from .f and .c files.
#
.f.o:
		$(FC) $(FFLAGS) $(IFLAGS) -c $*.f

.c.o:
		$(CC) $(CFLAGS) $(IFLAGS) -c $*.c
#
# C. If you just type 'make' this is equivalent to 'make all' where 'all'
#    means 'getorb lodr'. 'make test' runs lots of tests.
#
all:		$(EXECUTABLES)
test:		getorb_test lodr_test
		@echo "****************************************"
		@echo "*** Congratulations: all tests passed **"
		@echo "****************************************"
#
# D. Create the program sysdep to evaluate system dependencies and
#    the files created by it
#
sysdep:		underscore.c sysdep.f
		$(CC) -c underscore.c
		$(FC) $(FFLAGS) sysdep.f -o $@ underscore.o $(LFLAGS)
		rm -f underscore.o sysdep.o
sysdep.h sysdep.mk:	sysdep
		./sysdep
#
# E. Some OBJS need includes.
#
getorb.o:	math.inc
geocen.o:	GRS80.inc
fastio.o:	sysdep.h
#
# F. Now to test the software, type 'make tests'.
# 
getorb_test:	getorb
		@echo "****************************************"
		@echo "*** Testing getorb on ODR.378        ***"
		@echo "****************************************"
		./getorb t=950222025000,950222025100,1 data > getorb.out
		$(DOS2UNIX) getorb.out
		diff getorb.out data/getorb.out
		@echo "****************************************"
		@echo "*** Files are identical. Test passed ***"
		@echo "****************************************"
		rm -f getorb.out

lodr_test:	lodr
		@echo "****************************************"
		@echo "*** Testing lodr on ODR.378          ***"
		@echo "****************************************"
		./lodr t=950222.1,950222030000 data/ODR.378 > lodr.out
		$(DOS2UNIX) lodr.out
		diff lodr.out data/lodr.out
		@echo "****************************************"
		@echo "*** Files are identical. Test passed ***"
		@echo "****************************************"
		rm -f lodr.out
#
# G. To install the executables in the 'bin' directory.
#
install:	$(EXECUTABLES)
		chmod 755 $(EXECUTABLES)
		mv $(EXECUTABLES) $(BIN_DIR)
#
# H. If that is finished, cut the crap and type 'make clean'.
#    Or really get rid off all that was needed to make the executables
#    and type 'make spotless'.
#
clean:
		rm -f *.o sysdep sysdep.h sysdep.mk core
spotless:	clean
		rm -f `cat .cvsignore`
#
# I. We may like to make an archive
#
tar:		getorb.tar.gz
getorb.tar.gz:
		(cd .. ; tar -cvhf - $(SOURCES) | gzip > $@)
#
# J. Finally how to make the getorb library and executables
#
getorb.a:	$(OBJS)
		rm -f $@
		$(AR) $(ARFLAGS) $@ $(OBJS)
		$(RANLIB) $@
getorb:		getorb.a getorb_main.o
		$(FC) getorb_main.o getorb.a -o $@ $(LFLAGS)
lodr:		getorb.a lodr.o
		$(FC) lodr.o getorb.a -o $@ $(LFLAGS)
mdate:		mdate_main.o
		$(FC) mdate_main.o getorb.a -o $@ $(LFLAGS)
#
# K. Extra: convdate. Since it requires interaction with C, we
#    do it separately.
#
convdate:	convdate.o strf1970.o strf1985.o
		$(FC) convdate.o strf1970.o strf1985.o getorb.a -o $@ $(LFLAGS)
strf1970.o:	sysdep.h
