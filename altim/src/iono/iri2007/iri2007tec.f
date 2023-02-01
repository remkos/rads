**IRI2007TEC -- Compute TEC from IRI2007 model at two altitudes
*+
      SUBROUTINE IRI2007TEC (JMAG, UTC, LAT, LON, H1, H2, TEC1, TEC2)
      INTEGER*4 JMAG
      REAL*8	UTC, LAT, LON, H1, H2, TEC1, TEC2

* This routine computes the total electron content based on the
* IRI2007 model up to two altitudes (H1 and H2).
* Both altitudes should be above the F2 layer. If not, zeros
* are returned.
* Other input parameters are the geographical or geomagnetic
* latitude and longitude and the time in seconds since 1 Jan 1985.
*
* The code uses routines from the original IRI2007 code.
*
* Input arguments:
*   JMAG : 0 = geographic, 1 = geomagnetic coordinates
*   UTC  : UTC time in seconds since 1 Jan 1985
*   LAT  : Latitude (deg)
*   LON  : Longitude (deg)
*   H1   : Lower altitude (m)
*   H2   : Higher altitude (m)
*
* Output arguments:
*   TEC1 : Total electron content until H1 (TEC units)
*   TEC2 : Total electron content until H2 (TEC units)
*
* (1 TEC unit = 10^16/m^2)
*-
* $Log: iri2007tec.f,v $
* Revision 1.4  2013/03/26 23:17:30  rads
* - Add numbers to warning about layer heights.
*
* Revision 1.3  2009/01/22 18:51:32  rads
* - Guard against NaN input
*
* Revision 1.2  2008/01/21 16:21:42  rads
* - Made code more robust in several places
*
* Revision 1.1  2008/01/18 19:43:19  rads
* - Initial revision of this code
*
* Copyright (c) Remko Scharroo, Altimetrics LLC
*-----------------------------------------------------------------------
      real*4	hr(7),hmf2,tec,xlat,xlon,hour
* Highest accuracy
*     real*4	step(6)/1.0,0.5,1.0, 1.0, 1.0, 1.0/
* Recommended larger stepsize (within 1 permille)
      real*4	step(6)/2.0,1.0,2.5,10.0,30.0,30.0/
* Even larger stepsize (within 2 permille)
*     real*4	step(6)/4.0,2.0,5.0,20.0,50.0,50.0/
      integer*4	i,mjd,yy,mm,dd
      real*8	xmjd
      common	/block1/ hmf2
      save	step

* Guard against invalid inputs

      if (isnan(utc) .or. isnan(lat) .or. isnan(lon) .or.
     |	isnan(h1) .or. isnan(h2)) then
         tec1 = 0d0
	 tec1 = tec1/tec1
	 tec2 = tec1
         return
      endif

* Initialize

      tec1=0d0
      tec2=0d0
      xmjd=utc/86400d0+46066d0
      mjd=int(xmjd)
      call mjd2ymd(mjd,yy,mm,dd)
      hour=(xmjd-mjd)*24d0+25d0
      xlat=lat
      xlon=lon

* Call IRISUB to set up profile

      call irisub(jmag,xlat,xlon,yy,mm*100+dd,hour)

* Set various altitudes. Make sure both h1 and h2 are far above F layer.

      hr(1) = 100.
      hr(2) = hmf2-10.
      hr(3) = hmf2+10.
      hr(4) = hmf2+150.
      hr(5) = hmf2+250.
      hr(6) = h1*1e-3
      hr(7) = h2*1e-3
      if (hr(6).lt.hr(5)) then
         write (0,*) "iri2007tec: h1 should be above F2 layer:", hr(6), hr(5)
	 return
      else if (hr(7).lt.hr(6)) then
         write (0,*) "iri2007tec: h2 should be larger than h1:", hr(7), hr(6)
	 return
      endif

* Do the numerical integration until H1 and H2

      do i=1,6
         call iri2007tec1(hr(i),hr(i+1),step(i),tec)
	 tec2=tec2+tec
      enddo
      tec1=tec2-tec
      return
      end

      subroutine iri2007tec1(h1,h2,step,tec)
      implicit none
      real*4	h1,h2,step,tec,xe,dh,h
      integer*4 i,n
      if (h2.le.h1) then
         tec=0.0
         return
      endif
      n=max(1,nint((h2-h1)/step))
      tec=0.5*(xe(h1)+xe(h2))
      dh=(h2-h1)/n
      h=h1
      do i=1,n-1
         h=h+dh
         tec=tec+xe(h)
      enddo
      tec=tec*dh*1e-13
      return
      end
