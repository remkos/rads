        program ascii2grd

* This program convert the ASCII distributed files for
* the GOT Tide Model to DEOS grid files in order to save
* disk space and time in reading the files.
*-
*  1-Oct-2001 - Remko Scharroo
* 10-Jun-2003 - Create grids of different scale
* 25-Oct-2004 - Made able to deal with airtide grids
*-----------------------------------------------------------------------
      implicit none
      integer*4 nxmax,nymax,iunit,i,j,k,l,lnblnk
      parameter (nxmax=1441,nymax=721)
      integer*4 nx,ny,dn,mx,fdim,fdre,openf,iarg,iargc
      real*8 xmin,ymin,xmax,ymax
      real*8 dx, dy, mask,pi,rad,re,im,hmin,hmax,fact,scale
      real*8 wra(nxmax,nymax), wrg(nxmax,nymax)
      integer*2 iim(nxmax), ire(nxmax), i2min
      parameter (i2min=-32768)
      character*20 format
      character head*160,filenm*80,unit*2
      logical ltlend

      if (iargc().eq.0) then
	 write (*,500)
	 goto 9999
      endif
500   format('ascii2grd_got -- Convert GOT Ascii to DEOS grid file'//
     |'syntax: ascii2grd_got <ascii_files>')
      iunit=10
      pi=4*atan(1d0)
      rad=pi/180

* Read input parameters

      do iarg=1,iargc()
	 hmin=1d30; hmax=-1d30

	 call getarg(iarg,filenm)
         l=lnblnk(filenm)
         print *,' ASCII file : ', filenm(1:l)
         open (iunit,file=filenm,status='old',err=100)

* Read ascii file

550	 format (a)
	 read(iunit,550,err=101) head(1:80)
	 read(iunit,550,err=101) head(81:160)
         read(iunit,*,err=101) ny,nx
         read(iunit,*,err=101) ymin,ymax
         read(iunit,*,err=101) xmin,xmax
         read(iunit,*,err=101) mask
	 read(iunit,550,err=101) format
	 dx=(xmax-xmin)/(nx-1)
	 dy=(ymax-ymin)/(ny-1)

         do j=1,ny
            read(iunit,format) (wra(i,j),i=1,nx)
         enddo
	 do j=1,7
	    read(iunit,*,err=101)
	 enddo
         do j=1,ny
            read(iunit,format) (wrg(i,j),i=1,nx)
         enddo

         close(iunit) 

* Write grid file

	 l=index(head,' ')-1
         print *,' Binary files : ', head(:l)//'_{im,re}.grd'
	fdim=openf(head(:l)//'_im.grd','w')
	 if (fdim.lt.0) goto 102
	fdre=openf(head(:l)//'_re.grd','w')
	 if (fdre.lt.0) goto 102

         mx=nx/2*2+1
600      format ('@GRID',2e14.7,2i7,1x,4f8.4)
	 if (head(:8).eq.'M2 Ocean') then
	    fact=1.5d-2
	 else
	    fact=1d-2
	 endif
	 if (index(head,'(cm)').gt.0) then
	    scale=fact*1d-2
	    unit='cm'
	 else if (index(head,'(mm)').gt.0) then
	    scale=fact*1d-3
	    unit='mm'
	 else if (index(head,'(microbar)').gt.0) then
	    fact=1d-1
	    scale=fact*1d-3
	    unit='microbar'
	 else
	    stop 'can not figure out unit of input grid values'
	 endif
         write (head,600) 0d0,scale,mx,ny,0.0,360.0,ymin,ymax

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
	 print *,' Conversion ready. Min/Max value (',unit,',-) = ',
     |		hmin,hmax,nint(hmin/fact),nint(hmax/fact)
          
         goto 110
 100     print *,' Problem opening ASCII file ... please check'
         goto 110
 101     print *,' Problem opening header file ... please check'
         goto 110
 102     print *,' Problem writing BINARY file ... please check'
 110     continue
 
      enddo
9999  end
