**CMAPRD -- Read colormap
*+
      INTEGER FUNCTION CMAPRD (NAME, IC0, IC1, RGB)
      CHARACTER*(*) NAME
      INTEGER IC0, IC1
      REAL RGB(3,0:*)
*
* Read colormap from file specified by FILENM. If this file does not
* exist, the file is read from a directory indicated by the environment
* variable CMAP_DIR.
* If the file is not found at either place, IC0 and IC1 will be 0.
*
* Arguments:
*  NAME    (input): Name of the colormap file.
*  IC0    (output): Begin of the colormap for surface plotting.
*  IC1    (output): End of the colormap for surface plotting.
*  RGB    (output): RGB values for each color index. E.g. RGB(2,N) will
*                   contain the GREEN component of index N.
*                   All components range from 0.0-1.0.
*  CMAPRD (output): error code: 0 = no error, 1 = error opening file,
*		    2 = format error.
*-
*  7-Apr-1993 - Created [Remko Scharroo]
* 29-Apr-1994 - New error codes
*-----------------------------------------------------------------------
      logical exist
      character*80 filenm,cmapdir
      character*4 spec
      integer*2 jc,jc0,jc1,jr,jg,jb
      integer l,ic

      cmaprd=0
      filenm=name
      ic0=0
      ic1=0
      inquire (file=filenm,exist=exist)
      if (.not.exist) then
	 call getenv('CMAP_DIR',cmapdir)
	 l=index(cmapdir,' ')-1
	 filenm=cmapdir(:l)//'/'//name
      endif
      open (92,file=filenm,status='old',form='unformatted',
     .   access='direct',recl=8,err=1310)
      read (92,rec=1) spec,jc0,jc1
      if (spec.ne.'@COL') goto 1320
      ic0=jc0
      ic1=jc1
      l=1
   10 l=l+1
      read (92,rec=l,iostat=ic) jc,jr,jg,jb
      if (ic.ne.0) return
      ic=jc
      rgb(1,ic)=real(jr)/32640
      rgb(2,ic)=real(jg)/32640
      rgb(3,ic)=real(jb)/32640
      goto 10

 1310 cmaprd=1
      return
 1320 cmaprd=2
      return
      end
