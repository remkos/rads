**PGPXCONT -- Plot contours in pixel array
*+
      SUBROUTINE PGPXCONT (NX, NY, MX, X, I, C, C0, C1, DC, IC)
      INTEGER NX, NY, MX, I(MX,*), IC
      REAL    X(MX,*), C(MX,*), C0, C1, DC

      integer kx,ky

* Convert value array X to contour-level array C.

      do ky=1,ny
         do kx=1,nx
	    if (x(kx,ky).gt.1e20) then
	       c(kx,ky)=-1
	    else
	       c(kx,ky)=(max(c0,min(x(kx,ky),c1))-c0)/dc
	    endif
         enddo
      enddo

* Determine whether contour passes between two adjacent pixels

      do ky=1,ny-1
         do kx=1,nx-1
	    call pgpxcnt0(c(kx,ky),c(kx+1,ky),c(kx,ky+1),i(kx,ky),ic)
	 enddo
      enddo

      end

      subroutine pgpxcnt0(c1,c2,c3,i1,ic)
      real c1,c2,c3,dc
      integer i1,ic

      if (c1.lt.0) return
      dc=sqrt((c1-c2)**2+(c1-c3)**2)

* If contours are close, do not plot them all.

      if (dc.gt.1.2) then
	 if (mod(int(c1),4).ne.0) return
      else if (dc.gt.0.6) then
	 if (mod(int(c1),2).ne.0) return
      endif
      if (c2.ge.0 .and. int(c1).ne.int(c2)) i1=ic
      if (c3.ge.0 .and. int(c1).ne.int(c3)) i1=ic
      end
