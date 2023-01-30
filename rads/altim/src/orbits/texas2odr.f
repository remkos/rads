      program texas2odr

* texas2odr: Convert UT/CSR POE format to DEOS ODR format
* 
* syntax: texas2odr poefilename odrfilename
*
* The POE format looks like:
*
* ....+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8....+....9....+....0....+....1....+....2....+....3
* X01 UTC     60.0 19920413140200.000000 19920518140200.000000 BODY-FIXED   6378136.300 298.257000 UT/CSR ERS-1 TEG-3 PRELIMINARY     
*  2448726.000000000E+00  0.84722222222222E-01                        0.21387953416249E+07 -0.11579982948715E+07  0.67284673558010E+07
*   0.39095459704462E+04 -0.60347544498949E+04 -0.22766137982632E+04  0.70235573854134E+02  0.33156775663164E+03  0.79534296923929E+06
*  2448726.000000000E+00  0.85416666666666E-01                        0.23674867264941E+07 -0.15186111764783E+07  0.65788051591236E+07
*   0.37111083288578E+04 -0.59814942907160E+04 -0.27105091771612E+04  0.66974946489920E+02  0.32732208385572E+03  0.79480363954961E+06
* ....+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8....+....9....+....0....+....1....+....2....+....3
*
* The ODR format is specified in ODR.man
*
* 04-Mar-1997 - Created by Remko Scharroo

      implicit none
      real*8    ae0,ae,f,epoch,dh,lat,lon,hgt
      real*8	pos(3),vel(3),jd,frac
      parameter (ae0=6378137.0d0)
      integer*4 i,nrec,beg,arc,rep,ver,itime,iargc,unit
      character arg*80,spec*3,line*132

      if (iargc().ne.2) then
	 write (6,600)
	 goto 9999
      endif
600   format ('texas2odr: convert UT/CSR POE format to',
     |' DEOS ODR format'//
     |'syntax: texas2odr poefilename odrfilename')

* Open POE file

      call getarg(1,arg)
      if (arg.eq.'-') then
         unit=5
      else
         unit=9
         open (unit,file=arg,form='formatted',status='old')
      endif

* Open ODR

      call getarg(2,arg)
      open (10,file=arg,form='unformatted',status='new',
     |	recl=16,access='direct')
      nrec=0
      rep=35000
      ver=600
      beg=0
      arc=0
      i=index(arg,'ODR.')
      if (i.gt.0) read (arg(i+4:i+6),*) arc

* Read POE header

550   format (a)
      read (unit,501) spec,epoch,ae,f
501   format (a3,16x,f19.6,36x,f11.3,f11.6)
      if (spec.ne.'X01') stop "texas2odr: input not POE file"

* Check parameters

      if (f.ne.298.257d0) stop "wrong flattening"
      dh=ae-ae0

* Convert records
* Check for possible second header: if found, ignore first data record

10    read (unit,550,end=100) line
      if (line(:3).eq.'X01') then
         read (unit,550)
         read (unit,550)
         read (unit,550) line
      endif
      read (line,510,end=100) jd,frac,pos
      read (unit,511,end=100) vel,lat,lon,hgt
510   format (d22.9,d22.14,22x,3d22.14)
511   format (6d22.14)
      itime=nint((jd-2446066.5d0+frac)*86400)
      if (mod(itime,60).eq.0) then
	 nrec=nrec+1
	 if (lon.gt.180) lon=lon-360
	 write (10,rec=nrec+2)
     |		itime,nint(lat*1d7),nint(lon*1d7),nint((hgt+dh)*1d3)
      endif
      goto 10

* Close ODR file

100   continue
      write (10,rec=1) 'xODRERS-2   ',beg
      write (10,rec=2) rep,arc,nrec,ver
      write (6,*) arc,nrec
      close (unit)
      close (10)

9999  end
