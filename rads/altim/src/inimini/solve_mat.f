      subroutine solve_mat
      implicit none
      include "data.inc"
      include "stat.inc"
      include "satcat.inc"
      real*4	w
      integer*4	it,ip,im,icol,vpnt,mpnt

* Add constraints to normal matrix

      do it=1,ngood
         w=max(orberr(satel(it)),min(sigma(it),
     |          maxerr(satel(it))))
*	 w=t_nres(it)/w**2/npar
	 w=1/w**2
         do icol=1,npar
            ip=vpnt(it,icol)
            im=mpnt(it,it,icol,icol)
            matrix(im)=matrix(im)+w
         enddo
      enddo
      
* Solve normal equations using factorisation

      call spbtrf('U',bndlng,bndwth-1,matrix,bndwth,it)
      call spbtrs('U',bndlng,bndwth-1,1,matrix,bndwth,vector,bndlng,it)
      end
