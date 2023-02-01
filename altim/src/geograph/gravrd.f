*+GRAVRD - Read gravity field (GEODYN II style)
**
      subroutine gravrd(fact,filenm)
      real*4 fact
      character*(*) filenm
*-
* 18-Mar-1993 - Created (Remko Scharroo)
* 12-Apr-1998 - Add initialisation and normal equation support
*  9-Aug-2001 - Move initialisation of h() to d_lmp
*----------------------------------------------------------------------
      include "geograph.inc"
      integer l,m,k,ltot,mtot,ip,jp,jcs,jdeg,jord,lnblnk,ind(npar)
      real*8 c,s
      real*4 d(npar),scale(npar)
      character com*6,line*80
      logical donorm,first/.true./
      save first

* Initialize if not already done

      if (first) then
         gm=398600.4415d9
         ae=6378136.3d0
         g(1)=0
         do ip=2,npar
            g(ip)=g(ip-1)+(ip-1)
         enddo
         do ip=1,npar
	    cs(ip)=0
	    sig(ip)=0
	    ics(ip)=0
	    ideg(ip)=0
	    iord(ip)=0
	    rvect(ip)=0
         enddo
         do ip=1,npar*(npar+1)/2
            covar(ip)=0
         enddo
	 ipmax=0
         if (test) write (*,*) 'gravinit: initialisation done'
      endif

* If input file is covariance or normal matrix:
* Store the covariance/normal matrix

      l=lnblnk(filenm)
      docovar=(filenm(l-3:l).eq.'.COV'.or.filenm(l-3:l).eq.'.cov')
      donorm =(filenm(l-3:l).eq.'.NOR'.or.filenm(l-3:l).eq.'.nor')

      if (docovar .or. donorm) then
         open (31,file=filenm,status='old',form='unformatted')
         read (31) m
         do ip=1,m
	    read (31) jcs,jdeg,jord,d(1),d(2),d(3)

* Find a slot for this parameter
* (1) in the initialisation phase (first) use the slot k directly
* (2) otherwise, search for the right slot

	    if (jdeg.gt.lmax) then
	       k=0
	    else if (first) then
	       ipmax=ipmax+1
               if (ipmax.gt.npar) then
	          write (0,*) 'too many parameters (ipmax>npar):',
     |		ipmax,npar
	          stop
	       endif
	       k=ip
	    else
	       do k=1,ipmax
	       if (
     |jcs.eq.ics(k).and.jdeg.eq.ideg(k).and.jord.eq.iord(k)) goto 2
               enddo
	       if (test) write (*,*)
     |		'gravrd: no match found:',jcs,jdeg,jord
	       k=0
2              continue
	    endif
	    if (k.gt.0) then
	       cs(k)=cs(k)+fact*solve(k)
	       scale(k)=d(3)
	       if (first) then
	          solve(k)=d(1)*d(3)
	          ics(k)=jcs
	          ideg(k)=jdeg
	          iord(k)=jord
	       endif
	          sig(k)=d(2)*d(3)
*	    write (*,*) k,ics(k),ideg(k),iord(k),d(1),d(2)
	    endif
	    ind(ip)=k
         enddo
      endif
      if (docovar) then
         do ip=1,m
	    read (31) (d(jp),jp=1,ip)
	    k=ind(ip)
	    do jp=1,ip
	       l=ind(jp)
	       if (k.ne.l) d(jp)=d(jp)*2
	       if (k.eq.0 .or. l.eq.0) then
	       else if (k.lt.l) then
		  covar(g(l)+k)=covar(g(l)+k)+fact*d(jp)*scale(k)*scale(l)
	       else
		  covar(g(k)+l)=covar(g(k)+l)+fact*d(jp)*scale(l)*scale(k)
	       endif
	    enddo
         enddo
         goto 9000
      else if (donorm) then
         do ip=1,m
	    read (31) (d(jp),jp=1,ip)
	    k=ind(ip)
	    do jp=1,ip
	       l=ind(jp)
	       if (k.eq.0 .or. l.eq.0) then
	       else if (k.lt.l) then
		  covar(g(l)+k)=covar(g(l)+k)+fact*d(jp)/scale(k)/scale(l)
	       else
		  covar(g(k)+l)=covar(g(k)+l)+fact*d(jp)/scale(k)/scale(l)
	       endif
	    enddo
         enddo
	 read (31,end=9000) (d(ip),ip=1,m)
	 do ip=1,m
	    k=ind(ip)
	    if (k.ne.0) rvect(k)=rvect(k)+fact*d(ip)
	 enddo
         goto 9000
      endif

* Standard (ASCII) file
* Open file and read header

      open (31,file=filenm,status='old')
550   format (a)
610   format (a6,8x,2i3,4x,d20.8,d15.8,d13.1)
620   format (a6,2i2,2d15.3)
      rewind (31)
      read (31,550) line
      read (31,550) line
      read (line,610) com,ltot,mtot,gm,ae,flat
      ltot=min(ltot,ndeg)
      mtot=min(mtot,ndeg)
*     write (*,*) filenm,ltot,mtot

* Read data lines
* GCOEF = C and S, GCOEFC = C and sigma-C, GCOEFS = S and sigma-S
* SGCOEF = sigma-C and sigma-S

10    read (31,550,end=20) line
      com=line(1:6)
      if (com.eq.'SGCOEF') then
	 read (line,620) com,l,m,c,s
	 if (l.gt.lmax) goto 10
	 call addcoef(first,l,m,1,fact*c)
	 call addcoef(first,l,m,2,fact*s)
      else if (com.eq.'GCOEF ') then
         read (line,610) com,l,m,c,s
	 if (l.gt.lmax) goto 10
	 call addcoef(first,l,m,1,fact*c)
	 call addcoef(first,l,m,2,fact*s)
      else if (com.eq.'GCOEFC') then
         read (line,610) com,l,m,c
	 if (l.gt.lmax) goto 10
	 call addcoef(first,l,m,1,fact*c)
      else if (com.eq.'GCOEFS') then
         read (line,610) com,l,m,s
	 if (l.gt.lmax) goto 10
	 call addcoef(first,l,m,2,fact*s)
      endif
      if (test) write (*,*) "C",l,m,ip-1,ics(ip-1),cs(ip-1)
      if (test) write (*,*) "S",l,m,ip,ics(ip),cs(ip)
      goto 10

20    continue

9000  first=.false.
      close (31)
      end

      subroutine addcoef(first,l,m,i,c)
      logical first
      integer ip,l,m,i
      real*8  c
      include "geograph.inc"
      if (m.eq.0 .and. i.eq.2) return
      if (.not.first) then
         do ip=1,ipmax
	    if (ideg(ip).eq.l .and. iord(ip).eq.m .and. ics(ip).eq.i) then
	       cs(ip)=cs(ip)+c
	       return
	    endif
	 enddo
      endif
      ipmax=ipmax+1
      if (ipmax.gt.npar) stop 'too many parameters (ip>npar)'
      ideg(ipmax)=l
      iord(ipmax)=m
      ics(ipmax)=i
      cs(ipmax)=cs(ipmax)+c
      end
