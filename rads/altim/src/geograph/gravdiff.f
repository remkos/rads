      program gravdiff
      integer i
      character*80 model
      include "geograph.inc"
      lmax=70
      ipmax=npar
      call getarg(1,model)
      call gravrd(1.0,model)
      call gravrd(-1.0,'egm96')
      model="EGM96.70.NOR"
      call gravrd(0.0,model)
      do i=1,ipmax
         write (*,'(4i4,2d24.15)') i,ideg(i),iord(i),ics(i),
     |		-cs(i),-cs(i)/sig(i)
      enddo
      end
