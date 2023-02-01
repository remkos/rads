*
* create list of conspicious data from residual file.
*
* usage: badres gii-residual-file edit-sigma
*
      integer maxobs,maxpas,maxsta
      parameter(maxobs=200,maxpas=2500,maxsta=9999)
      real*8 times(maxobs,maxpas),resid(maxobs,maxpas),
     |sigma(maxobs,maxpas),deriv(maxobs,maxpas),elev(maxobs,maxpas)
      integer ipass(4,maxpas)
      character cpass(maxpas)*8
      real*8 t1985,delt,fact,resi,sigm
      parameter(t1985=(46066-30000)*86400d0)

      character*80 arg,filenm/' '/,date1*12,date2*12
      real*8 edit/3.5d0/
      integer ios,i,ip,ista,isat,nobs,type,ndel,iargc,it
      logical elim/.true./
c
      do i=1,iargc()
	 call getarg(i,arg)
	 if (arg(1:2).eq.'-h') then
	    write (0,610)
	    goto 9999
	 else if (arg(1:2).eq.'-e') then
	    elim=.true.
	 else if (arg(1:2).eq.'-t') then
	    elim=.false.
	 else if (arg(1:5).eq.'edit=') then
	    read (arg(6:),*) edit
	 else
	    filenm=arg
	 endif
      enddo

      if (filenm.eq.' ') filenm='fort.19'
      call resread(filenm,maxobs,maxpas,ipass,cpass,times,resid,sigma,
     |deriv,elev)
c
      ndel=0
      if (elim) write (6,603)
      do ip=1,maxpas
	 ista=ipass(1,ip)
	 isat=ipass(2,ip)
	 nobs=ipass(3,ip)
	 type=ipass(4,ip)
	 if (type.eq.40) then
	    delt=27d0
	    fact=1000d0
	 else
	    delt=0d0
	    fact=1d0
	 endif
	 do i=1,nobs
	    resi=resid(i,ip)*fact
	    sigm=sigma(i,ip)*fact
	    if (abs(resi).gt.sigm*edit) then
	       ndel=ndel+1
	       it=nint(times(i,ip)-t1985+delt)
	       if (elim) then
		  call strf1985(date1,'%y%m%d%H%M%S',it-1)
		  call strf1985(date2,'%y%m%d%H%M%S',it+1)
		  write (6,602) date1,date2,ista,isat,resi,sigm
	       else
		  call strf1985(date1,'%y%m%d%H%M%S',it)
		  call strf1985(it,date1)
		  write (6,600) date1,ista,isat,resi,sigm
	       endif
	    endif
	 enddo
      enddo
600   format (a12,i6,i9,2f12.3)
602   format (2(a12,2x),i4,i9,2f10.3)
603   format ('1')
610   format ('badres: list bad residuals from residual file'//
     |'usage: badres [options] [residual-file]'//
     |'where [residual-file] is the name of the GII residual file',
     |' (default: fort.19)'/
     |'and   [options] are'/
     |'      edit=edit : mark all residuals edit*sigma as bad',
     |' (default: 3.5)'/
     |'      -e        : produce input for ELIMREC (default)'/
     |'      -t        : produce list of trouble')
c
9999  end
