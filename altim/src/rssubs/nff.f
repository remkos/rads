**NFF -- Check status of NetCDF return code
*+
      SUBROUTINE NFF(IOS)
      INTEGER*4 IOS

* This routine checks the return value of any NetCDF call and verifies
* if an error had resulted. If so, the program will terminate.
*
* Example:
*     call nff(nf_open(filenm,nf_nowrite,ncid))
*-
* $Log: nff.f,v $
* Revision 1.2  2006/06/05 17:11:06  rads
* - Initial version
*
*-----------------------------------------------------------------------
      include "netcdf.inc"
      if (ios.ne.nf_noerr) call fin("NetCDF error: "//nf_strerror(ios))
      end
