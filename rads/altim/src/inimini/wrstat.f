      subroutine wrstat
      implicit none
      include "stat.inc"
      include "data.inc"
      include "init.inc"
      include "satcat.inc"
      logical first/.true./
      integer n,ix,iy,ip
      character*4 sat1,sat2
      save first
      real edit,w,erms

* On first call print statistics header

      if (first) then
	 first=.false.
	 write (6,1001)
	 do ix=1,maxsat
	    call shortsat(ix,sat1)
	    if (s_ncor(ix).ne.0) write (6,1010) sat1
	 enddo
	 if (inimode.eq.1) then
	    do ix=1,maxsat
	       do iy=1,ix
	          ip=h(ix)+iy
	          call shortsat(iy,sat1)
	          call shortsat(ix,sat2)
	          if (s_nres(ip).ne.0) write (6,1020) sat1,'x',sat2
	       enddo
	    enddo
	 else
	    do ix=1,maxsat
	       ip=h(ix)+ix
	       call shortsat(ix,sat1)
	       if (s_nres(ip).ne.0) write (6,1020) sat1,'-','ref='
	    enddo
	 endif
	 write (6,550) '|== total ==|=tracks=|'
      endif

* Write iteration number

      write (6,1000) iter

* Write mean and rms of the orbit correction per satellite

      adjust=0
      do ix=1,maxsat
	 n=s_ncor(ix)
	 if (n.ne.0) then
	    s_mcor(ix)=s_mcor(ix)/n
	    s_rcor(ix)=s_rcor(ix)/n-s_mcor(ix)**2
	    if (s_rcor(ix).gt.0) s_rcor(ix)=sqrt(s_rcor(ix))
	    write (6,1030) nint(s_mcor(ix)*1000),
     |		nint(s_rcor(ix)*1000),
     |		nint(s_adjust(ix)*1000)
	    adjust=max(adjust,s_adjust(ix))
	 endif
      enddo

* Write number and rms of residuals per satellite combination

      nrestot=0
      resrms=0e0
      do ip=1,maxcmb
	 n=s_nres(ip)
	 if (n.ne.0) write (6,1040) n,
     |		nint(sqrt(s_rres(ip)/n)*1000)
	 nrestot=nrestot+n
	 resrms=resrms+s_rres(ip)
      enddo
      resrms=sqrt(resrms/nrestot)

* Determine track sigma and good tracks
*
* Tracks are marked 'bad' when
* - Number of used residuals is less or equal to the number of estimated
*   parameters
* - Residual rms per track > edttrk times the overall residual rms
*   (edttrk best value = 2?)
*
* Track sigma equals the residual rms per track (divided by sqrt(2)
* for xovers)

      erms=0
      ip=1
      do ix=1,ngood
	 ip=ip+t_nres(ix)
	 erms=erms+t_eres(ix)
      enddo
      erms=sqrt(erms/ip)

      ip=0
      do ix=1,ngood
         w=sqrt(t_rres(ix)/t_nres(ix))

         good(ix)=.true.
	 edit=w/resrms
	 edit=100.-100.*t_eres(ix)/t_nrestot(ix)
	 if (edit.gt.edttrk
     |		.and. iter.gt.0) then
            good(ix)=.false.
	    write (99,'("Rej 1:",2i10,f10.6)') iter,itrinv(ix),edit
         else if (t_nres(ix).le.2*npar) then
	    good(ix)=.false.
	    write (99,'("Rej 2:",3i10)') iter,itrinv(ix),t_nres(ix)
	 endif
*	 do icol=1,npar
*	    it=vpnt(ix,icol)
*	    if (param(it).gt.maxerr(satel(ix))) then
*	    write (99,'("Rej 3:",2i10,f10.6)') iter,itrinv(ix),param(it)
*	    good(ix)=.false.
*	    endif
*	 enddo
         if (good(ix)) then
	    ip=ip+1
	 else
	    w=1d3
	 endif

	 if (inimode.eq.1) w=w/sqrt(2.)
         w=min(w,maxerr(satel(ix)))
*        w=max(orberr(satel(ix)),min(w,maxerr(satel(ix))))
         sigma(ix)=w
      enddo

      ntrks=ip
      write (6,1040) nrestot,nint(resrms*1000)
      write (6,1050) ntrks
      if (iter.eq.0) resrms0=resrms

550   format (a)
1000  format (i2,1x,$)
1001  format (/'Processing statistics:'/
     .'- per sat : orbit adjustments in mm (mean, rms,',
     .' max in the last iteration)'/
     .'- per sat-sat combination : residuals (number and rms in mm)'//
     .'   ',$)
1010  format ('|====',a,'====',$)
1020  format ('|=',a,a,a,'=',$)
1030  format (i5,i4,i4,$)
1040  format (i8,i4,$)
1050  format (i8)

      end

      subroutine shortsat(ix,sat)
      integer ix,l,lnblnk
      character sat*4,satname*8
      include "satcat.inc"

      satname=name(ix)
      sat=satname(1:4)
      l=lnblnk(satname)
      if (satname(l:l).ge.'0' .and. satname(l:l).le.'9')
     |	 sat(4:4)=satname(l:l)
      end
