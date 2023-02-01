        program ascii2grd

* This program convert the ASCII distributed files for
* FES Tide Model to DEOS grid files in order to save
* disk space and time in reading the files.
*-
* First version : Jean-Marc Molines 15/11/1994
* Revised       : Fabien Lefevre 09/08/2000
* DEOS grids    : Remko Scharroo 01/10/2001
* 10-Jun-2003 - Changed scale to 0.1 mm for all grids except M2 (0.15 mm)
* 19-Nov-2004 - Increased grid size to allow for FES2004
*-----------------------------------------------------------------------
      implicit none
      integer*4 nxmax,nymax,iunit,i,j,k,l,lnblnk,iargc,iarg
      parameter (nxmax=2881, nymax=1441)
      integer*4 nx,ny,dn,mx,fdim,fdre,openf,nwave
      real*8 xmin,ymin,xmax,ymax
      real*8 dx,dy,mask,pi,rad,re,im,hmin,hmax,fact
      real*8 wra(nxmax,nymax),wrg(nxmax,nymax)
      integer*2 iim(nxmax), ire(nxmax), i2min
      parameter (i2min=-32768,nwave=14)
      character*126 file_asc,file_bin
      character*80 head
      logical	ltlend

      if (iargc().eq.0) then
	 write (*,500)
	 goto 9999
      endif
500   format('ascii2grd_fes -- Convert FES Ascii to DEOS grid file'//
     |'syntax: ascii2grd_fes <ascii_files>')
      iunit=10
      pi=4*atan(1d0)
      rad=pi/180

* Read input parameters

      call system ('mkdir -p ../load ../ocean')

      do iarg=1,iargc()
	 hmin=1d30; hmax=-1d30

	 call getarg(iarg,file_asc)
	 l=index(file_asc,'_')-1
	 if (index(file_asc,'drfes').gt.0) then
	    file_bin='../load/'//file_asc(:l)
	 else
	    file_bin='../ocean/'//file_asc(:l)
	 endif
	 if (index(file_asc,'M2_fes').gt.0) then
	    fact=1.5d-2
	 else
	    fact=1d-2
	 endif

         l=lnblnk(file_asc)
         print *,' Reading ASCII file ... be patient ...'
         print *,' File : ', file_asc(1:l)
         open (iunit,file=file_asc,status='old',err=100)

* Read ascii file

         read(iunit,*,err=101) xmin,xmax
         read(iunit,*,err=101) ymin,ymax
         read(iunit,*,err=101) dx,dy
         read(iunit,*,err=101) nx,ny
         read(iunit,*,err=101) mask,mask

         do j=1,ny
            do k=1,nx,30
              read(iunit,*) (wra(i,j),i=k,MIN(nx,k+29))
              read(iunit,*) (wrg(i,j),i=k,MIN(nx,k+29))
            enddo
         enddo

         close(iunit) 

* Write grid file

         print *,' Writing grid files ... will not be so long ...'
         l=lnblnk(file_bin)
	 fdre=openf(file_bin(:l)//'_re.grd','w')
	 if (fdre.lt.0) goto 102
	 fdim=openf(file_bin(:l)//'_im.grd','w')
	 if (fdim.lt.0) goto 102

         mx=nx/2*2+1
600      format ('@GRID',2e14.7,2i7,1x,4f8.4)
         write (head,600) 0d0,fact/1d2,mx,ny,0d0,360d0,ymin,ymax

         call writef(fdim,80,head)
         call writef(fdre,80,head)

* Shift grid to range from 0 to 360 degrees

         if (xmin.eq.-180) then
	    dn=nx/2
	 else if (xmin.eq.0) then
	    dn=0
	 else
	    stop "grid should start at -180 or 0"
	 endif
         do j=1,ny
           do i=1,nx
             if(i.le.dn) then
               k=i+dn
             else
               k=i-dn
             endif
             if((wra(i,j).eq.mask).or.(wrg(i,j).eq.mask)) then
               iim(k)=i2min
               ire(k)=i2min
             else
               re=wra(i,j)*cos(wrg(i,j)*rad)
               im=wra(i,j)*sin(wrg(i,j)*rad)
	       hmin=min(hmin,im,re)
	       hmax=max(hmax,im,re)
               iim(k)=nint(im/fact)
	       ire(k)=nint(re/fact)
             endif
           enddo
           if (ltlend()) then
	      call i2swap(nx,iim)
	      call i2swap(nx,ire)
           endif
	   iim(mx)=iim(1)
	   ire(mx)=ire(1)
	   call writef(fdim,2*mx,iim)
	   call writef(fdre,2*mx,ire)
         enddo

         call closef(fdim)
         call closef(fdre)
	 print *,' Conversion ready. Min/Max value (cm,-) = ',
     |		hmin,hmax,nint(hmin/fact),nint(hmax/fact)
          
         goto 110
 100     print *,' Problem opening ASCII file ... please check'
         goto 110
 101     print *,' Problem opening header file ... please check'
         goto 110
 102     print *,' Problem writing BINARY file ... please check'
 110     continue
	 fact=1d-2
 
      enddo
9999  end
