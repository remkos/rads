      subroutine whatindx(idx)
* Determine index 
      integer kx,ky,idx
      include 'xrgb.inc'

      kx=x/xblksize+1
      ky=(y-.25)/yblksize+1
      idx=(9-ky)*16+(kx-1)
      if (kx.ge.1 .and. kx.le.16 .and. ky.ge.1 .and. ky.le.9) return
      do idx=-12,-1,1
	 if (x.ge.bx0(-idx) .and. x.le.bx1(-idx) .and.
     .		y.ge.by0(-idx) .and. y.le.by1(-idx)) return
      enddo
      idx=-9
      end
