**IRI95TEC -- Compute TEC from IRI95 model at two altitudes
*+
      SUBROUTINE IRI95TEC (JMAG, UTC, LAT, LON, H1, H2, TEC1, TEC2)
      implicit none
      INTEGER*4 JMAG
      REAL*8	UTC, LAT, LON, H1, H2, TEC1, TEC2

* This routine computes the total electron content based on the
* IRI95 model up to two altitudes (H1 and H2).
* Both altitudes should be above the F2 layer. If not, zeros
* are returned.
* Other input parameters are the geographical or geomagnetic
* latitude and longitude and the time in seconds since 1 Jan 1985.
*
* The code uses routines from the original IRI95 code.
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
* 28-Jun-2002 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      real*4	hr(7),hmf2,xnmf2,hmf1,tec
* Highest accuracy
*     real*4	step(6)/1.0,0.5,1.0, 1.0, 1.0, 1.0/
* Recommended larger stepsize (within 1 permille)
*     real*4	step(6)/2.0,1.0,2.5,10.0,30.0,30.0/
* Even larger stepsize (within 2 permille)
      real*4	step(6)/4.0,2.0,5.0,20.0,50.0,50.0/
      integer*4	i,mjd,yymmdd,mdate
      real*8	xmjd
      common	/block1/ hmf2,xnmf2,hmf1
      save	step

* Initialize

      tec1=0d0
      tec2=0d0
      xmjd=utc/86400d0+46066d0
      mjd=int(xmjd)
      yymmdd=mdate(3,mjd)

* Call IRISUB to set up profile

      call irisub(jmag,real(lat),real(lon),yymmdd/10000,
     |	mod(yymmdd,10000),real((xmjd-mjd)*24+25))

* Set various altitudes. Make sure both h1 and h2 are far above F layer.

      hr(1) = 100.
      hr(2) = hmf2-10.
      hr(3) = hmf2+10.
      hr(4) = hmf2+150.
      hr(5) = hmf2+250.
      hr(6) = h1/1e3
      hr(7) = h2/1e3
      if (hr(6).lt.hr(5)) then
         write (*,*) "iri95tec: h1 should be above F2 layer"
	 return
      else if (hr(7).lt.hr(6)) then
         write (*,*) "iri95tec: h2 should be larger than h1"
	 return
      endif

* Do the numerical integration until H1 and H2

      do i=1,6
         call iri95tec1(hr(i),hr(i+1),step(i),tec)
	 tec2=tec2+tec
      enddo
      tec1=tec2-tec
      return
      end

      subroutine iri95tec1(h1,h2,step,tec)
      implicit none
      real*4	h1,h2,step,tec,xe,dh,h
      integer*4 i,n
      n=nint((h2-h1)/step)
      if (n.le.0) then
         tec=0.0
         return
      endif
      dh=(h2-h1)/n
      tec=(xe(h1)+xe(h2))/2
      h=h1
      do i=1,n-1
         h=h+dh
         tec=tec+xe(h)
      enddo
      tec=tec*dh/1e13
      return
      end
