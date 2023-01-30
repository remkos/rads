      program stapos

* Program to convert station positions:
* - to and from XYZ and Lat-Lon-Height
* - add or subtract eccentricities in XYZ or North-East-Up
*-
* 10-May-2001 - Created from genstapos
*-----------------------------------------------------------------------
      implicit	none
      integer*4	msta
      parameter (msta=1000)
      integer*2 pnt(9999)
      real*8	s_pos(6,msta),s_sig(6,msta),xyz1(3),xyz2(3)
      integer*4	s_id(msta),s_plate(msta),noise,s_wave(msta)
      real*8	s_ecc(3,msta),ecc(3),dxyz(3),epoch,dum,rad,
     |		r1,lat1,lon1,hgt1,r2,lat2,lon2,hgt2
      integer*4	i,j,nsta,statinfo,mjdref,option,lnblnk,id/9999/,ios,
     |		way/2/
      character text*80,xyzfile*80,s_name(msta)*32,s_code(msta)*4
      logical	neu/.true./

* Initialise

      do i=1,9999
         pnt(i)=0
      enddo
      rad=atan(1d0)/45
      xyzfile='/user/remko/ers/setups/stations.xyz'
      epoch=1d30

* Read arguments

      call getarg(1,text)
      if (text(:2).eq.'-h') then
         write (*,1300)
	 goto 9999
      else if (text.ne.' ') then
         xyzfile=text
      endif
1300  format ('stapos -- manipulate station coordinates'//
     |'syntax: stapos [coordinatefile]'//
     |'where:'/
     |' coordinatefile : file with station coordinates',
     |' (def: /user/remko/setups/stations.xyz)')

* Read the XYZ coordinates and velocities from the default file

      call loadcoord(xyzfile,mjdref,nsta,s_id,s_pos,s_sig)

* Set pointers (only set those that are zero)
* Reduce velocities from m/year to m/s.
* Default adjustment to zero.
* Load system.data info

      do i=1,nsta
         if (pnt(s_id(i)).eq.0) pnt(s_id(i))=i
	 s_plate(i)=9876
	 j=statinfo(epoch,s_id(i),s_name(i),dum,s_code(i),
     |		noise,s_wave(i),s_plate(i),s_ecc(1,i))
	 if (j.ne.0) write (*,1301) s_id(i)
      enddo

* Add dummy station

      nsta=nsta+1
      s_id(nsta)=9999
      pnt(9999)=nsta
      s_name(nsta)='Generic station'
      s_code(nsta)=' '
      do i=1,3
         s_ecc(i,nsta)=0
      enddo

* Start interactive part of program

550   format (a)
551   format (a,' -> ',$)

* Ask for station coordinates

25    write (*,551)
     |'Enter coordinates (XYZ or LatLonHgt) or station id'
      read (*,550) text
      if (lnblnk(text).le.4) then
	 read (text,*) id
	 i=pnt(id)
	 if (i.eq.0) then 
	    write (*,1302) id
	    goto 25
	 endif
	 do j=1,3
	    xyz1(j)=s_pos(j,i)
	    ecc(j)=s_ecc(j,i)
	 enddo
      else
         read (text,*) xyz1
	 ecc(1)=0;ecc(2)=0;ecc(3)=0
      endif

* Coordinates:
* xyz1 : monument coordinates
* xyz2 : optical centre coordinates
* dxyz : eccentricity

* Check if XYZ or LatLonHgt

      call xyzorllh(xyz1,r1,lat1,lon1,hgt1)

* Convert eccentricity N-E-Up to XYZ (if required)

40    continue
      if (neu) then
        dxyz(1)=-ecc(1)*sin(lat1)*cos(lon1)-ecc(2)*sin(lon1)
     |		+ecc(3)*cos(lat1)*cos(lon1)
        dxyz(2)=-ecc(1)*sin(lat1)*sin(lon1)+ecc(2)*cos(lon1)
     |		+ecc(3)*cos(lat1)*sin(lon1)
        dxyz(3)=+ecc(1)*cos(lat1)+ecc(3)*sin(lat1)
      endif

* Way=1: compute monument from optical centre and eccentricity
* Way=2: compute optical centre from monument and eccentricity
* Way=3: compute eccentricity from monument and optical centre

      if (way.eq.1) then
         do j=1,3
	    xyz1(j)=xyz2(j)-dxyz(j)
         enddo
         call xyzgeo(xyz1,r1,lat1,lon1,hgt1)
      else if (way.eq.2) then
         do j=1,3
	    xyz2(j)=xyz1(j)+dxyz(j)
         enddo
         call xyzgeo(xyz2,r2,lat2,lon2,hgt2)
      else
         do j=1,3
	    dxyz(j)=xyz2(j)-xyz1(j)
         enddo
	 neu=.false.
      endif

* Convert eccentricity XYZ to N-E-Up (if required)

      if (.not.neu) then
         ecc(1)=(lat2-lat1)*r1
         ecc(2)=(lon2-lon1)*r1*cos(lat1)
         ecc(3)=hgt2-hgt1
      endif

* Print station info

      write (*,600) id,s_code(i),s_name(i)
      write (*,610) 'Monument coordinates',xyz1,lat1/rad,lon1/rad,hgt1
      write (*,611) 'Eccentricity',dxyz,ecc
      write (*,610) 'Optical centre coordinates',xyz2,
     |lat2/rad,lon2/rad,hgt2
600   format (/
     |'Station nr and code  : ',i4,1x,a/
     |'Station location     : ',a)
610   format (/a/
     |'  X   Y   Z      [m] :',3f14.4/
     |' Lat Lon Hgt [deg,m] :',2f14.7,f14.4)
611   format (/a/
     |' dX  dY  dZ      [m] :',3f14.4/
     |' dN  dE  dU      [m] :',3f14.4)

* Ask for eccentricity

552   format (/'Enter eccentricity (',a,
     |') [',3f10.4,'] -> ',$)
555   format (/'Enter ',a,' coordinates [',3f14.4,'] -> ',$)
554   format (/'Enter station id [',i4,'] -> ',$)
553   format (/'What do you want to do?'/
     |' 0. Quit'/
     |' 1. Change eccentricity and recompute monument coordinates'/
     |' 2. Change eccentricity and recompute optical',
     |' centre coordinates'/
     |' 3. Change monument coordinates and recompute optical',
     |' centre coordinates'/
     |' 4. Change monument coordinates and recompute eccentricity'/
     |' 5. Change optical centre coordinates and recompute',
     |' monument coordinates'/
     |' 6. Change optical centre coordinates and recompute',
     |' eccentricity'/
     |' 7. Start over'/
     |'Enter option [7] -> ',$)
560   format (/'Eccentricity is given:'/
     |' 1. by N-E-Up'/' 2. by X-Y-Z'/' 3. by station id'/
     |'Enter option [3] -> ',$)
60    write (*,553)
      read (*,550,end=9990) text
      option=7
      read (text,*,iostat=ios) option
      neu=.true.
      if (option.le.0 .or. option.gt.7) then
	 goto 9999
      else if (option.le.2) then
         write (*,560)
	 read (*,550) text
	 i=3
	 read (text,*,iostat=ios) i
	 if (i.eq.1) then
            write (*,552) 'N-E-Up',ecc
            read (*,550) text
            read (text,*,iostat=ios) ecc
         else if (i.eq.2) then
            write (*,552) 'X-Y-Z',dxyz
            read (*,550) text
            read (text,*,iostat=ios) dxyz
	    neu=.false.
         else
	    write (*,554) id
	    read (*,550) text
            read (text,*,iostat=ios) id
            i=pnt(id)
            if (i.eq.0) then 
	       write (*,1302) id
	       goto 60
            endif
            do j=1,3
	       ecc(j)=s_ecc(j,i)
            enddo
         endif
	 way=option
      else if (option.le.4) then
         write (0,555) 'monument',xyz1
         read (*,550) text
         read (text,*,iostat=ios) xyz1
	 call xyzorllh(xyz1,r1,lat1,lon1,hgt1)
	 way=2
	 if (option.eq.4) way=3
      else if (option.le.6) then
         write (0,555) 'optical centre',xyz2
         read (*,550) text
         read (text,*,iostat=ios) xyz2
	 call xyzorllh(xyz2,r2,lat2,lon2,hgt2)
	 way=1
	 if (option.eq.6) way=3
      else
         goto 25
      endif
      goto 40

1301  format ('WARNING: Station ',i4.4,' not found in system.data. ',
     |'Station ignored. Please edit system.data')
1302  format ('WARNING: Station ',i4.4,' not found in station ',
     |'coordinate file.')

9990  write (*,550)
9999  end

      subroutine xyzorllh(xyz,r,lat,lon,hgt)
      real*8 xyz(3),r,lat,lon,hgt,rad
      rad=atan(1d0)/45
      if (xyz(1)**2+xyz(2)**2+xyz(3)**2.lt.6d6**2) then
	 lat=xyz(1)*rad
	 lon=xyz(2)*rad
	 hgt=xyz(3)
	 call geoxyz(lat,lon,hgt,xyz,r)
      else
	 call xyzgeo(xyz,r,lat,lon,hgt)
      endif
      end
