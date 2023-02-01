#
# Makefile for gdr2mar models
#
# Specify Fortran and C compiler commands. In most "make"
# distributions these variables are already set to the
# appropriate values.
# To set them to the right values, change the next lines and
# remove the # character at the beginning of the line.
#
#FC     = g77
#CC	= gcc
#
# Specify the flags to be used with the Fortan and C compilers.
# You will need to figure out whether you need the
# -DUNDERSCOREBEFORE, -DUNDERSCOREAFTER or -DCAPITALS option
# in order to use the C routines of fastio.c in the gdr2mar.f
# program. The options for some architectures are given below.
# Note the special provisions for HP Unix compilers!
#
#FFLAGS = -O
CFLAGS	= -O -DUNDERSCOREAFTER	# Linux and IRIS
#CFLAGS	= -O			# AIX and HPUX
#FFLAGS	= -O +U77		# HPUX
#LFLAGS	= -lU77			# HPUX
#
# NO CHANGES should be needed in the next lines
#
all:	gdr2mar

gdr2mar:	gdr2mar.o fastio.o lnblnk.o
	$(FC) -o $@ gdr2mar.o fastio.o lnblnk.o $(LFLAGS)

clean:
	rm -f gdr2mar gdr2mar.o fastio.o lnblnk.o

tar:
	( cd .. ; tar -cvzhf gdr2mar.tar.gz gdr2mar )
