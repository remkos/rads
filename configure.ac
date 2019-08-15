# -----------------------------------------------------------------------
#  Copyright (c) 2011-2019  Remko Scharroo
#  See LICENSE.TXT file for copying and redistribution conditions.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
# -----------------------------------------------------------------------
#
# Run GNU autoconf to produce a configure script from this input.
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# Define and show package version
# -----------------------------------------------------------------
AC_INIT([RADS], [4.3.6], [remko.scharroo@eumetsat.int])
AC_PREREQ([2.67])
AC_MSG_NOTICE([AC_PACKAGE_NAME, version AC_PACKAGE_VERSION])

# -----------------------------------------------------------------
# Check for this file before proceeding
# -----------------------------------------------------------------
AC_CONFIG_SRCDIR(src/rads.f90)

# -----------------------------------------------------------------
# Set some defaults
# -----------------------------------------------------------------
AC_PREFIX_DEFAULT(`pwd`)
AC_CANONICAL_SYSTEM

# -----------------------------------------------------------------
# Specify user options
# -----------------------------------------------------------------
AC_ARG_WITH([nf-config],  [  --with-nf-config=dir       Full path to nf-config (or nc-config) [[auto]]])
AC_ARG_WITH([netcdf-inc], [  --with-netcdf-inc=dir      Specify location of the NetCDF include file [[auto]]])
AC_ARG_WITH([netcdf-lib], [  --with-netcdf-lib=dir:dir  Specify location of the NetCDF library.  Use colon-separator if Fortran library is in separate directory [[auto]]])
AC_ARG_ENABLE([debug],    [  --enable-debug             Compile with debugging options [[off]]])
AC_ARG_ENABLE([devel],    [  --enable-devel             Compile developer tools [[off]]])

# -----------------------------------------------------------------
# Process user options
# -----------------------------------------------------------------
AC_MSG_NOTICE([checking user options])
configure_options="`$as_echo "$ac_configure_args" | sed "s/^ //; s/\'//g"`" #'
AC_SUBST(CONFIGURE_OPTIONS, $configure_options)

# -----------------------------------------------------------------
# Switch for debugging
# -----------------------------------------------------------------
AC_MSG_CHECKING([whether debugging is desired])
if ! test "x$enable_debug" = "xyes" ; then
	enable_debug=no
	ac_cv_prog_cc_g=no
	ac_cv_prog_fc_g=no
fi
AC_MSG_RESULT([$enable_debug])

# -----------------------------------------------------------------
# Switch for developers only
# -----------------------------------------------------------------
AC_MSG_CHECKING([whether developer tools are desired])
if ! test "x$enable_devel" = "xyes" ; then
	enable_devel=no
fi
AC_MSG_RESULT([$enable_devel])

# -----------------------------------------------------------------
# Default settings come from nf-config
# -----------------------------------------------------------------
AC_PATH_PROGS([nf_config], [nf-config] [nc-config], [$with_nf_config])

# -----------------------------------------------------------------
# Find Fortran 90 compilers.
# -----------------------------------------------------------------
# If F90 is set, and FC isn't, then set FC to whatever F90 was set to.
# Otherwise set it to what NetCDF was compiled with
AC_MSG_CHECKING([for Fortran 90 compiler])
if test "${F90+set}" = "set"; then
	F90=$F90
elif test "${FC+set}" = "set"; then
	F90=$FC
elif ! test "x$nf_config" = "x"; then
	F90=`$nf_config --fc`
	F90=`basename $F90`
fi
AC_MSG_RESULT([$F90])

# On solaris, f95 comes before gfortran.
case "$host" in
*solaris*) F90G="f95 f90 gfortran" ;;
	*) F90G="gfortran f95 f90" ;;
esac
F90S="xlf90 fort xlf95 ifort ifc efc pgf95 lf95 $F90G xlf90 pgf90 epcf90"

# Find a fortran 90 compiler, if requested.
AC_MSG_NOTICE([finding Fortran 90 compiler])
AC_PROG_FC($F90 $F90S)
F90=$FC
if test "${FC+set}" = set && test "x$FC" = x; then
	AC_MSG_ERROR([Cannot find F90 compiler. Cannot continue.])
fi

# -----------------------------------------------------------------
# Learn more about the Fortran 90 compiler
# -----------------------------------------------------------------
AC_MSG_NOTICE([setting up Fortran 90])
# Find out various additions needed to FCFLAGS to make the compiler
# run with .f90, free-form, 132-byte length, or module inclusion.
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM
AC_FC_LINE_LENGTH(132)
AC_FC_MODULE_FLAG
AC_FC_FUNC("FC_name",FC_name)
AC_LANG_POP(Fortran)
AC_SUBST(MOD_FLAG, ${ac_cv_fc_module_flag})
AC_SUBST(FC_name)

# Check if the fortran 90 compiler can compile a simple test file.
AC_MSG_CHECKING([if $FC can handle Fortran 90])
cat <<EOF >conftest.f90
program foo
real, parameter :: pi = 3.141592654
end program
EOF
doit='$FC $FCFLAGS conftest.f90 -o conftest.exe'
if AC_TRY_EVAL(doit); then
	AC_MSG_RESULT([yes])
else
	AC_MSG_ERROR([Cannot use $F90 to compile Fortran 90 code. Cannot continue.])
fi
rm -rf conftest.*

# -----------------------------------------------------------------
# Additional flags for GNU compilers (for developers only)
# -----------------------------------------------------------------
AC_MSG_CHECKING([for additional compiler flags (for developers)])
opt=
if test "x$enable_devel" = "xyes"; then
	if test "x$ac_cv_f77_compiler_gnu" = "xyes" ; then
		opt=" -Wall -Wimplicit"
	fi
	# Version 4.* uses -Wtabs, Version 5.* (and later?) uses -Wno-tabs
	if test "x$FC" = "xgfortran" ; then
		gfortran --version | grep -q " 4."
		if test $? == 1 ; then
			opt=" -Wall -Wno-tabs -fimplicit-none"
		else
			opt=" -Wall -Wtabs -fimplicit-none"
		fi
	elif test "x$FC" = "xifort" ; then
		opt=" -warn"
	fi
fi
FCFLAGS=$FCFLAGS$opt
AC_MSG_RESULT([$opt])

if test "x$enable_devel" = "xno"; then
	ALTIM=
	DEVEL_DIRS=
	DEVEL_OBJS=
	DEVEL_FLAGS=
	INSTALL_C=
else
	DEVEL_DIRS="devel"
	DEVEL_OBJS='${ALTIM}/src/rssubs/tpj_subs.o ${ALTIM}/src/rssubs/solar_subs.o'
	DEVEL_FLAGS=" -I../devel $MOD_FLAG"'${ALTIM}/include'
	INSTALL_C=" -C"
fi
AC_SUBST(ALTIM)
AC_SUBST(DEVEL_DIRS)
AC_SUBST(DEVEL_OBJS)
AC_SUBST(DEVEL_FLAGS)

# -----------------------------------------------------------------
# Allow user to set --with-netcdf-inc and/or --with-netcdf-lib
# -----------------------------------------------------------------
AC_MSG_CHECKING([for NetCDF includes])
if ! test "x$with_netcdf_inc" = "x" ; then
	NETCDF_INC=$MOD_FLAG$with_netcdf_inc
elif test -x "$nf_config" ; then
	NETCDF_INC=$MOD_FLAG`$nf_config --includedir`
else
	NETCDF_INC=$MOD_FLAG/usr/local/inc
fi
AC_MSG_RESULT([$NETCDF_INC])
AC_SUBST(NETCDF_INC)

AC_MSG_CHECKING([for NetCDF libraries])
if ! test "x$with_netcdf_lib" = "x" ; then
	NETCDF_LIB="-L"`echo $with_netcdf_lib | sed "s=:= -L="`" -lnetcdf"
	doit="find "`echo $with_netcdf_lib | sed "s=:= ="`" -name libnetcdff*"
	if test `$doit|wc -l` -ne 0; then
		NETCDF_LIB="${NETCDF_LIB}f -lnetcdf"
	fi
elif test -x "$nf_config" ; then
	NETCDF_LIB=`$nf_config --flibs`
else
	NETCDF_LIB=-L/usr/local/lib -lnetcdf -lnetcdff
fi
AC_MSG_RESULT([$NETCDF_LIB])
AC_SUBST(NETCDF_LIB)

# -----------------------------------------------------------------
# Try to compile simple program with NetCDF
# -----------------------------------------------------------------
AC_MSG_CHECKING([if $FC can compile a simple NetCDF program])
cat <<EOF > conftest.f90
program foo
use typesizes
use netcdf
integer :: i = 0, j
write (*,"(a)") nf90_inq_libvers()
if (i == 1) j = nf90_def_var_deflate (0, 0, 0, 0, 0)
end program
EOF
doit="$FC $FCFLAGS conftest.f90 -o conftest.exe $NETCDF_INC $NETCDF_LIB"
if AC_TRY_EVAL(doit); then
	AC_MSG_RESULT([yes])
	NETCDF_VERS=`./conftest.exe`
	AC_MSG_RESULT([Succeeded to compile Fortran 90 code with NetCDF version $NETCDF_VERS])
	rm -f ./conftest.exe
else
	AC_MSG_WARN([no])
	echo $doit
	echo "========================================================"
	echo "Could not compile simple Fortran 90 program with NetCDF."
	echo "One or more of the next problems occurred:"
	echo "- The directory/directories containing libnetcdf/libnetcdff"
	echo "  library/libraries were not/wrongly identified. Use the"
	echo "  --with-netcdf-lib=DIR[:DIR] option to specify the proper"
	echo "  directory/directories."
	echo "- The directory containing netcdf.mod was not/wrongly"
	echo "  identified. Use the --with-netcdf-inc=DIR option to"
	echo "  specify the proper directory."
	echo "- The NetCDF library pointed to was not compiled with"
	echo "  Fortran 90 support. Recompile the NetCDF library with the"
	echo "  configuration flag --enable-f90"
	echo "- Your compiler $FC does not properly compile Fortran 90 code."
	echo "  Find another compiler and specify it on the command line"
	echo "  using FC=<compilername>."
	echo "- Your NetCDF Fortran library is too old. Use version 4.3"
	echo "  or later."
	echo "========================================================"
	exit 1
fi

# -----------------------------------------------------------------
# Allow user to set --with-lapack-lib
# -----------------------------------------------------------------
AC_MSG_CHECKING([for lapack options (for developers)])
LAPACK_LIB=
if test "x$with_devel" = "xyes" ; then
	case $build_os in
	*darwin*) LAPACK_LIB="-Wl,-framework,Accelerate" ;;
	*) LAPACK_LIB="-llapack -lblas" ;;
	esac
fi
AC_MSG_RESULT([$LAPACK_LIB])
AC_SUBST(LAPACK_LIB)

# -----------------------------------------------------------------
# Do we have big or little endian
# -----------------------------------------------------------------
AC_C_BIGENDIAN
if test "x$ac_cv_c_bigendian" = "xyes"; then
	LITTLE_ENDIAN="false"
else
	LITTLE_ENDIAN="true"
fi
AC_SUBST(LITTLE_ENDIAN)

# -----------------------------------------------------------------
# Find other utilities
# -----------------------------------------------------------------
AC_MSG_NOTICE([finding other utilities])
AC_PROG_LN_S
AC_PROG_RANLIB
AC_PROG_INSTALL

# -----------------------------------------------------------------
# Add -C option to install (developers only)
# -----------------------------------------------------------------
INSTALL_DATA="${INSTALL_DATA}${INSTALL_C}"
INSTALL_PROGRAM="${INSTALL_PROGRAM}${INSTALL_C}"
INSTALL_SCRIPT="${INSTALL_SCRIPT}${INSTALL_C}"
INSTALL_DIR='${INSTALL} -d'
AC_SUBST(INSTALL_DIR)

# -----------------------------------------------------------------
# Substitute in .in files
# -----------------------------------------------------------------
AC_OUTPUT(config.mk src/rads-config doc/manuals/rads_config.tex)
