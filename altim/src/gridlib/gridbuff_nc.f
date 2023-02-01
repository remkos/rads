**GRIDBUFF_NC -- Load a grid into memory (NetCDF format)
*+
      FUNCTION GRIDBUFF_NC (FILENM, POINTER)
      INTEGER*4 GRIDBUFF_NC, POINTER
      CHARACTER FILENM*(*)
*-
* 26-Jul-2006 - Avoiding use of %val
* 21-Aug-2005 - Created by Remko Scharroo from GRIDBUFF
*-----------------------------------------------------------------------
      integer*4 i,j,l,lnblnk,nbytes,mallocf,pntr2,nb,noff,
     |		x_id,y_id,z_id,nvars,dims(2),ncid,ndims,memloc
      real*8	dummy(2)
      logical	nfs
      include "gridbuff.inc"
      include "nan.inc"
      include "netcdf.inc"

* Implicit definition of nfs()

      nfs(i)=(i.ne.nf_noerr)

* Open grid file

      l=lnblnk(filenm)
      if (nfs(nf_open(filenm,nf_nowrite,ncid))) goto 1310

* Interpret header information
* Look for first 2-dimensional (z) variable and determine

      if (nfs(nf_inq_nvars(ncid,nvars))) goto 1320
      do i=1,nvars
         if (nfs(nf_inq_varndims(ncid,i,ndims))) goto 1320
	 if (ndims.eq.2) then
	    z_id=i
	    goto 50
	 endif
      enddo
      write (0,1300) 'Could not find 2-D variable in',filenm(:l)
      goto 1330
50    continue

* Get the data type

      if (nfs(nf_inq_vartype(ncid,z_id,ntype))) goto 1320
      if (ntype.eq.1) then
         nb=1
      else if (ntype.eq.3) then
         nb=2
      else if (ntype.eq.4) then
         nb=4
      else if (ntype.eq.5) then
         nb=4
      else if (ntype.eq.6) then
         nb=8
      else
         write (0,1300) 'Unknown data type in',filenm(:l)
	 goto 1330
      endif

* Get the ids of the x and y variables

      if (nfs(nf_inq_vardimid(ncid,z_id,dims))) goto 1320
      x_id=-1
      y_id=-1
      do i=1,nvars
	 if (nfs(nf_inq_varndims(ncid,i,ndims))) goto 1320
	 if (ndims.eq.1) then
	    if (nfs(nf_inq_vardimid(ncid,i,j))) goto 1320
	    if (j.eq.dims(1)) x_id=i
	    if (j.eq.dims(2)) y_id=i
	 endif
      enddo
      if (x_id.lt.0 .or. y_id.lt.0) then
	 write (0,1300) 'Could not find the x or y variables in',filenm(:l)
         goto 1330
      endif
      if (nfs(nf_inq_dimlen(ncid,dims(1),nx))) goto 1320
      if (nfs(nf_inq_dimlen(ncid,dims(2),ny))) goto 1320

* Get z-range, -scale and -offset and missing value

      if (nfs(nf_get_att_double(ncid,z_id,"scale_factor",dz))) dz=1d0
      if (nfs(nf_get_att_double(ncid,z_id,"add_offset"  ,z0))) z0=0d0

      if (nfs(nf_get_att_double(ncid,z_id,"actual_range",dummy)) .and.
     |    nfs(nf_get_att_double(ncid,z_id,"valid_range" ,dummy))) then
         dummy(1)=-1d40
	 dummy(2)=1d40
      endif
      zmin=dummy(1)
      zmax=dummy(2)

      if (nfs(nf_get_att_double(ncid,z_id,"_FillValue"   ,znan)) .and.
     |	  nfs(nf_get_att_double(ncid,z_id,"missing_value",znan)))znan=nan

* Get x- and y-ranges

      if (nfs(nf_get_att_int(ncid,nf_global,"node_offset",noff))) noff=0
      if (nfs(nf_get_att_double(ncid,x_id,"actual_range",dummy)) .and.
     |    nfs(nf_get_att_double(ncid,x_id,"valid_range" ,dummy))) then
         if (nfs(nf_get_var1_double(ncid,x_id,   0,dummy(1)))) goto 1320
         if (nfs(nf_get_var1_double(ncid,x_id,nx-1,dummy(2)))) goto 1320
	 noff=0
      endif
      dx=(dummy(2)-dummy(1))/(nx+noff-1)
      xmin=dummy(1)+noff*dx/2
      xmax=dummy(2)-noff*dx/2

      if (nfs(nf_get_att_double(ncid,y_id,"actual_range",dummy)) .and.
     |    nfs(nf_get_att_double(ncid,y_id,"valid_range" ,dummy))) then
         if (nfs(nf_get_var1_double(ncid,y_id,   0,dummy(1)))) goto 1320
         if (nfs(nf_get_var1_double(ncid,y_id,ny-1,dummy(2)))) goto 1320
	 noff=0
      endif
      dy=(dummy(2)-dummy(1))/(ny+noff-1)
      ymin=dummy(1)+noff*dy/2
      ymax=dummy(2)-noff*dy/2

* Determine dimensions and allocate memory

      nbytes=nx*ny*nb		! Number of bytes of grid data
      nbuf=(nbytes+mhead+7)/8*8	! Round up to multiple of 8 bytes
      if (mallocf(nbuf,pointer).gt.0) goto 1330
      pntr2=pointer+mhead-memloc(tmp_b)

* Store grid information in the data buffer header

      call memput(pointer,mhead,head)

* Load entire data block.

      if (ntype.eq.1) then
         if (nfs(nf_get_var_int1(ncid,z_id,tmp_b(pntr2)))) goto 1340
      else if (ntype.eq.3) then
         if (nfs(nf_get_var_int2(ncid,z_id,tmp_b(pntr2)))) goto 1340
      else if (ntype.eq.4) then
         if (nfs(nf_get_var_int(ncid,z_id,tmp_b(pntr2)))) goto 1340
      else if (ntype.eq.5) then
         if (nfs(nf_get_var_real(ncid,z_id,tmp_b(pntr2)))) goto 1340
      else
         if (nfs(nf_get_var_double(ncid,z_id,tmp_b(pntr2)))) goto 1340
      endif
      GRIDBUFF_NC=0
      goto 9999

* Error exits

1300  format ('GRIDBUFF_NC: ',a,1x,a)
1310  GRIDBUFF_NC=1
      goto 9998
1320  GRIDBUFF_NC=2
      goto 9998
1330  GRIDBUFF_NC=3
      goto 9998
1340  GRIDBUFF_NC=4

* Here normal termination

9998  POINTER=0
9999  i=nf_close(ncid)
      end
