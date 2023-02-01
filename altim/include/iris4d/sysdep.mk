F90	= f90
CPPFLAGS=#
FFLAGS	= -u -bytereclen -fullwarn -OPT:IEEE_comparisons=ON
FEXTSRC = -extend_source
RECUR	= -LANG:recursive=on
CPLUS	= -Xcpluscomm
MATHLIB	= -lcomplib.sgimath
SHARED_LD	= ld -shared -all
SHARED_LIB_LIBS =# -lX11 -lm -lc
SHARED_EXT	= so
