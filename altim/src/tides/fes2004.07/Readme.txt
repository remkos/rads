###############################################################
# Distributed software for the computation of tides 
# from the 2004 Finite Element Solution (FES2004)
###############################################################
#  This software have been tested on Linux platform
#  It is provided without any guarantees.
# 
#  For bug reports, please contact :
#  ---------------------------------
#  Fabien LEFEVRE 
# 
#  CLS
#  http://www.cls.fr
#  Direction Oc�anographie Spatiale
#  8-10, rue Herm�s - Parc Technologique du Canal
#  31526 Ramonville Saint-Agne cedex - France
#  Tel: +33 (0)5 61 39 37 45 Fax: +33 (0)5 61 39 37 82
#  e-mail: Fabien.Lefevre@cls.fr
# 
#  NOTE: This software is based on the former versions
#        developed by Jean-Marc MOLINES and Florent LYARD
###############################################################


-----------------
March, 22nd 2007
-----------------

This package is the fully revised version of the FES distribution with:
- new finite element meshes which were used to compute FES2004
- a new prediction algorithm debbuged and written in C
- a new F77 interface
- new loading effects computed from FES2004
- some bugs corrected on the S1 grid
- new K2 grids and associated loading effects

Please refer to file tesfes.c to get an example of the C routine to use

WE STRONGLY RECOMMAND TO USE THE C SOURCES
However for those who are familiar with Fortran, we added a F77/C interface.
Please refer to file fcc_tesfes.F to get an example of the fortran routine
to use.


---------------------
WHAT IS DISTRIBUTED ?
---------------------
In the ./ directory
- Readme.txt          : This file

In the doc/ directory you will find the following files:
- Readme.htm          : An html file which describes the Software
                        Development Kit (SDK) for the people who want
                        to add the FES package in their programs

In the src/ directory you will find the following files:
- Makefile            : Makefile to adapt to different platforms
- error.c             : Error management for FES algortihm
- fes.c               : Main package of subroutines for the algorithm
- fes-int.h           : Common declaration
- grid.c              : Package of subroutines to manage grids
- interp.c            : Package of subroutines to interpolate value
- prediction.c        : Package of subroutines for the prediction algorithm
- testfes.c           : Test program written in Cthe algorithm
- fcc_testfes.F       : A test program written in F77the algorithm
- fes.output.fes2004  : Output of testfes or fcc_testfes with FES2004
- ascii2bin.c         : Convert ascii FES files to binary fes files


In the include/ directory you will find the following files:
- fes.h               : Common declaration to include if the libfes.a is used

In the lib/ directory you will find the following files:
nothing : the libfes.a will be built in it 

In the data/ directory you will find the following files:
- 15 ascii files with the pure harmonic tides :
   - 2N2_fes2004.asc
   - K1_fes2004.asc
   - K2_fes2004.asc
   - M2_fes2004.asc
   - N2_fes2004.asc
   - O1_fes2004.asc
   - P1_fes2004.asc
   - Q1_fes2004.asc
   - S2_fes2004.asc
   - S1_fes2004.asc (thanks Richard Ray)
   - Mf_fes2004.asc
   - Mm_fes2004.asc
   - Mtm_fes2004.asc
   - Msqm_fes2004.asc
   - M4_fes2004.asc (thanks Thierry Letellier)
- 9 ascii files with the radial tide computed from FES2004 (thanks Olivier Francis) :
   - 2N2_drfes2004.asc
   - K1_drfes2004.asc
   - K2_drfes2004.asc
   - M2_drfes2004.asc
   - N2_drfes2004.asc
   - O1_drfes2004.asc
   - P1_drfes2004.asc
   - Q1_drfes2004.asc
   - S2_drfes2004.asc


--------------------------
HOW TO BUILD THE PACKAGE ?
--------------------------
Once the fes2004.tar file have been 'untared', you have to follow the 
following steps:

1/ Adapt the Makefile to the owner system (default is Linux)

2/ Select one (and only one) make command depending on the use of the
   algorithm:

   1) make CFLAGS="-O" ascii
      FES execution : load ascii files (long initialisation, few hard
      drive space)

   2) make CFLAGS="-O -DBINARY" binary
      the binary FES files will be generated once if needed
      FES execution : load binary files (need more hard drive space,
      lot of RAM on computer, but
      faster execution: useful for numerous tide computations)

   3) make CFLAGS="-O -DBINARY -DIO" binary
      the binary FES files will be generated once if needed
      FES execution : load binary files + allow direct access on the binary
      files (need more hard drive space, few RAM on computer,
      but very fast execution: useful for several tide computations) 

3/ The code is compiled; link and process a test:
   - testfes     : a C program to test the code.
   - fcc_testfes : a F77 program to test the code.

        
-------------------------
HOW TO USE THE SOFTWARE ?
-------------------------
Please follow the instructions written in:
    - testfes.c for the C interface
    - fcc_testfes.f for the F77 interface


---------------------
IN CASE OF PROBLEMS ?
---------------------
This software has been tested and ran successfully on Red Hat Linux release 7.1
(Seawolf) with gcc and pgf77 (be careful lot of problems occur when using g77!!!)

If you have problems, you can send me an email to Fabien.Lefevre@cls.fr

---------------------
BIBLIOGRAPHY
---------------------
Author: Lyard, F.; Lef�vre, F.; Letellier, T; Francis, O.
Year: 2006
Title: Modelling the global ocean tides: a modern insight from FES2004
Journal: Ocean Dynamics
Volume: 56
Issue: Numbers 5-6
Pages: 394-415
Date: Decembre 2006
Short Title: Modelling the global ocean tides: a modern insight from FES2004
Original Publication: Springer Berlin / Heidelberg



F. March, 22nd 2007
