      program ftest
      integer i
      integer*4 ndeg,narray
      parameter (ndeg=70,narray=((2*ndeg+3)*(ndeg+2)*(ndeg+1))/6-1)
      real*8    flmp(0:narray)
      integer*4 pnt(0:ndeg)
#ifdef DIFF
      real*8    dflmp(0:narray),ddflmp(0:narray)
#endif
      do i=1,10
#ifdef DIFF
         call f_lmpdd(ndeg,1d0,pnt,flmp,dflmp,ddflmp)
#else
         call f_lmp(ndeg,1d0,pnt,flmp)
#endif
      enddo
      do i=1,100
#ifdef DIFF
	 write (*,600) i,flmp(i),dflmp(i),ddflmp(i)
#else
	 write (*,600) i,flmp(i)
#endif
      enddo
600   format (i6,f16.11,f16.10,f16.8)
      end
