c
c     BLNEW is a subroutine that picks a new measurement
c     in the direction of RLATX. BLNEW = .TRUE. if error.
c
      function blnew(ibuf,asc,rdat,adat,irec,irec0,irec1,rlatx)
      integer*4 ibuf,irec,irec0,irec1,k
      real*8 rdat(5,3:4),rlatx
      integer*2 adat(3,3:4)
      logical blnew,asc

      blnew=.true.

      if (asc .eqv. (rdat(2,4).lt.rlatx)) then
         if (irec.ge.irec1) return
         do k=1,5
	    rdat(k,3)=rdat(k,4)
	 enddo
         irec=irec+1
         call blread(ibuf,irec+1,rdat(1,4),adat(1,4))
      else
         if (irec.le.irec0) return
         do k=1,5
	    rdat(k,4)=rdat(k,3)
	 enddo
         irec=irec-1
         call blread(ibuf,irec,rdat(1,3),adat(1,3))
      endif
      blnew=.false.
      end
