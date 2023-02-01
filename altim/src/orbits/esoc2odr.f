      program esoc2odr

* Convert ESOC orbit to ODR format
*
* usage: esoc2odr odrname < esocname
*
*-
*  2-Nov-1999 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4 yyyy,mm,dd,hh,mi,yyyymmdd,mdate,mjd,i,nrec,odr(4)
      real*8    ss,xyz(3),sec,r,lat,lon,height,rad
      character esocnm*80,odrnm*80,satel*8

* Read ESOC input and ODR output filenames. If omitted print usage information.

      call getarg (1,esocnm)
      call getarg (2,odrnm)
      if (odrnm.eq.' ') then
         write (0,600)
600      format ('esoc2odr - Convert ESOC orbit to ODR format'//
     |'usage: esoc2odr esocname odrname'//
     |'where:'/
     |'  esocname : name of the input ESOC file'/
     |'  odrname  : name of the output ODR file')
         goto 9999
      endif

* Initialize

      nrec=0
      if (index(esocnm,'.er1').gt.0) then
         satel='ERS-1'
      else if (index(esocnm,'.er2').gt.0) then
         satel='ERS-2'
      else
         write (0,"('Unknown satellite, assuming ERS-2')")
	 satel='ERS-2'
      endif
      rad=atan(1d0)/45d0

* Open the input ESOC file and output ODR file

      open (10,file=esocnm,status='old',form='formatted')
      open (20,file=odrnm,status='new',form='unformatted',
     |access='direct',recl=16)

* Read the ESOC file.
* Skip the header line and continue reading all lines until EOF.

      read (10,*)		! Skip header line
10    read (10,610,end=1000) yyyy,mm,dd,hh,mi,ss,xyz
610   format (i4,4(1x,i2),1x,f6.3,8x,4f15.6)

* Convert epoch to UTC seconds since 1.0 Jan 1985

      yyyymmdd=yyyy*10000+mm*100+dd
      mjd=mdate(2,yyyymmdd)
      sec=(mjd-46066)*86400d0+hh*3600d0+mi*60d0+ss

* Convert XYZ from kilometers to meters

      do i=1,3
         xyz(i)=xyz(i)*1d3
      enddo

* Convert Earth fixed XYZ to Geodetic coordinates

      call xyzgeo(xyz,r,lat,lon,height)

* Convert real values to integer values for ODR

      odr(1)=nint(sec)
      odr(2)=nint(lat/rad*1d7)
      odr(3)=nint(lon/rad*1d7)
      odr(4)=nint(height*1d3)

* Write ODR record

      nrec=nrec+1
      write (20,rec=nrec+2) odr
      goto 10

* End of conversion loop
* Write ODR headers and close file

1000  write (20,rec=1) 'xODR',satel,0
      write (20,rec=2) 35000,0,nrec,0
      close (20)

9999  end
