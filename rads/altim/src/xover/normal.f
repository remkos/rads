      subroutine normal(nrmpnt,hpnt,rpnt,rmsmax,sigma,rms,error)
      integer*4 nrmpnt
      real*8    hpnt(nrmpnt),rpnt(nrmpnt),rmsmax,sigma,rms
      logical   error

      integer   i,j,k,m,maxpnt,if
      parameter (maxpnt=61)
      real*8    a(maxpnt,4),ata(16)/16*0/,atr(4)/4*0/,
     |          v(4),h(4),hnrm,d
      integer   icall/0/,kp(4)
      save

      error=.true.
      if (icall.eq.0) then
	 icall=1
	 call matsy1(4,h)
	 do m=1,nrmpnt
	    do j=1,4
              a(m,j)=dble(m-nrmpnt/2-1)**(j-1)
	    enddo
	 enddo
	 k=0
	 do i=1,4
	    do j=1,i
	       k=k+1
	       do m=1,nrmpnt
                  ata(k)=ata(k)+a(m,j)*a(m,i)
	       enddo
	    enddo
	 enddo
	 if=1
	 call matcho(ata,h,v,kp,4,d,if)
       endif

       call mattmv(4,nrmpnt,maxpnt,a,hpnt,atr)
       call linch1(ata,h,atr,kp,4)
       call matmmv(nrmpnt,4,maxpnt,a,atr,rpnt)
       rms=0
       do m=1,nrmpnt
          rms=rms+(hpnt(m)-rpnt(m))**2
       enddo
       rms=sqrt(rms/nrmpnt)
       if (rms.gt.rmsmax) return
       do m=1,nrmpnt
          if (abs(hpnt(m)-rpnt(m)).gt.sigma*rms) return
       enddo
       error=.false.
       hnrm=rpnt(nrmpnt/2+1)
       end
