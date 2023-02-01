      subroutine fill_mat
      implicit none
      integer ip,im
      include "data.inc"
      logical first/.true./
      save first

* On first call: initialise memory

      if (first) then
	 first=.false.

* Determine
* - bandwidth of the matrix
* - size of the band matrix and vector
* If bndwth=0 use the largest possible value

         bndlng=ngood*npar
         if (bndwth.eq.0) bndwth=9999999
	 bndwth=min(bndlng,bndwth,maxsiz/bndlng)
	 bndwth=bndwth/npar*npar
	 if (inimode.eq.2) bndwth=npar
         write (6,1110) ngood,npar,bndwth,
     |		bndwth*bndlng*4,bndlng*4

      endif

* Initialise matrix

      do im=1,bndwth*bndlng
         matrix(im)=0e0
      enddo

* Initialise vector

      do ip=1,bndlng
         vector(ip)=0e0
      enddo

      return

 1110 format (/'Memory usage:'/
     |'Number of tracks in data file   :',i9/
     |'Number of parameters per track  :',i9/
     |'Half-width of the band matrix   :',i9/
     |'Size of the band matrix         :',i9,' bytes'/
     |'Size of the right-hand vector   :',i9,' bytes')
      end
