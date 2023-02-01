**FLATHER -- Give tidal height at given location according to Flather model
*+
      FUNCTION FLATHER (LAT, LON, MJD)
      REAL*8 FLATHER, LAT, LON, MJD
*
* This function returns the total Ocean Tidal elevation for a number of
* given tidal consituents, using the Flather model.
*
* Arguments:
*  LAT      (input): Geodetic latitude of the location in degrees.
*  LON      (input): East longitude of the location in degrees.
*  MJD      (input): Epoch.
*  FLATHER (output): Total Ocean Tidal elevation in metres.
*-
*  5-Jul-1993: New manual
*  4-Jul-1994: GRIDRD8 implemented
*  4-Oct-1994: All variables defined
*  5-Jan-1996: Declarations changed to data statement
*-----------------------------------------------------------------------
      integer*4 nxp,nx,kx,nyp,ny,ky,nc,k,icall,l
      real*8 x0,x1,x,dx,rx,y0,y1,y,dy,ry,sina,cosa,z0,z1
      character*80 filenm
      parameter (nxp=113,nyp=106,nc=8)
      real*8 h(nxp,nyp,nc),g(nxp,nyp,nc)
      real*8 phase(nc),amp(nc),tidht2
      character*2 compo(8)

      include 'flather.inc'

      data compo/'Q1','O1','P1','K1','N2','M2','S2','K2'/
      save

* Initialize on first call

      if (icall.eq.1) goto 100

      icall=1
      pi=4*atan(1d0)
      twopi=2*pi
      rad=pi/180
      nx=nxp
      ny=nyp
      filenm=directory
      l=index(filenm,' ')
      filenm(l:l)='/'
      do k=1,nc
	 filenm(l+1:l+2)=compo(k)
         filenm(l+3:l+4)='.H'
         call gridrd8(filenm,nx,ny,h(1,1,k),x0,x1,y0,y1,z0,z1)
	 if (nx.ne.nxp .or. ny.ne.nyp) goto 1320
         filenm(l+3:l+4)='.G'
         call gridrd8(filenm,nx,ny,g(1,1,k),x0,x1,y0,y1,z0,z1)
	 if (nx.ne.nxp .or. ny.ne.nyp) goto 1320
	 do kx=1,nxp
	    do ky=1,nyp
       		sina=h(kx,ky,k)*sin(g(kx,ky,k)*rad)
     		cosa=h(kx,ky,k)*cos(g(kx,ky,k)*rad)
		h(kx,ky,k)=sina
		g(kx,ky,k)=cosa
	    enddo
	 enddo
      enddo
      dx=(x1-x0)/(nx-1)
      dy=(y1-y0)/(ny-1)

* Check if requested point is within area boundaries

100   flather=1d35
      x=(lon-x0)/dx+1
      y=(lat-y0)/dy+1
      if (x.lt.1 .or. x.gt.nxp .or. y.lt.1 .or. y.gt.nyp) return

* Determine closest gridpoint to the bottom left

      kx=int(x)
      ky=int(y)

* Check if valid tidal height can be obtained.

      if (h(kx  ,ky,1)+h(kx  ,ky+1,1)+
     |    h(kx+1,ky,1)+h(kx+1,ky+1,1).gt.1d20) return
      rx=x-kx
      ry=y-ky

* Determine the Doodson phase for Greenwich.
* Add tidal sea height for several constituents to each corner.

      do k=1,nc
         sina=(1-rx)*((1-ry)*h(kx  ,ky,k)+ry*h(kx  ,ky+1,k))+
     |            rx*((1-ry)*h(kx+1,ky,k)+ry*h(kx+1,ky+1,k))
         cosa=(1-rx)*((1-ry)*g(kx  ,ky,k)+ry*g(kx  ,ky+1,k))+
     |            rx*((1-ry)*g(kx+1,ky,k)+ry*g(kx+1,ky+1,k))
         phase(k)=atan2(sina,cosa)/rad
         amp(k)=sqrt(sina**2+cosa**2)
      enddo
      flather=tidht2(mjd,lat,amp,phase)
      return

  550 format (a)
 1320 write (6,550) 'flather: error reading grid '//filenm
      end


      subroutine astro2 (ishunt,mjd,xlat)

      integer*4 ishunt
      real*8 mjd,xlat

*  IF ISHUNT IS 0, PARAMETERS FOR THE MAIN AND SATELLITE TIDAL
*  COMPONENTS ARE READ IN AND THE FREQUENCIES COMPUTED. IF ISHUNT IS 1,
*  THE ARRAYS F,V & U ARE SET UP/UPDATED FOR THE GIVEN TIME (PREFERABLY
*  THE 16TH OF THE MONTH) AND LATITUDE. IF ISHUNT IS 2 ONLY THE 
*  LATITUDE-DEPENDENT ARRAYS F & U ARE UPDATED.
*  MODIFIED FROM A PROGRAM BY MIKE FOREMAN, INSTITUTE OF OCEAN SCIENCES,
*  CANADA.

      integer*4 ncs,mjd1976
      parameter (ncs=10,mjd1976=42778)

      include 'flather.inc'

      integer*4   ii(ncomp),jj(ncomp),kk(ncomp),ll(ncomp),mm(ncomp),
     |          nn(ncomp),nj(ncomp),ldel(ncomp,ncs),mdel(ncomp,ncs),
     |          ndel(ncomp,ncs),ir(ncomp,ncs)
      real*8    semi(ncomp),ee(ncomp,ncs),ph(ncomp,ncs),
     |          uucos(ncomp,ncs),uusin(ncomp,ncs)
      integer*4 unit
      logical opened
      character*80 filenm
      integer*4 iv,k,j,l
      real*8 vdbl,uu,uudbl,enp,p,pp,dpp,dnp,rr,h,years,s,
     |  h0,s0,p0,enp0,pp0,dh,ds,dp,tau,dtau,slat,sums,sumc
      save 

*  FORMATS
1     format (5f13.10)
2     format (11x,1x,6i3,f5.2,i4)
3     format ((11x,3(3i3,f4.2,f7.4,1x,i1,1x)))

      if (ishunt .eq. 2) go to 666
      if (ishunt .eq. 1) go to 555

*  READ IN ASTRONOMICAL ARGUMENTS, DOODSON NUMBERS,PHASE CORRECTIONS
*  NO. OF SATELLITES FOR EACH COMPONENT AND CHANGES IN THE LAST THREE
*  DOODSON NUMBERS, PHASE CORRECTIONS, AMPLITUDE RATIO & FLAG FOR
*  LATITUDE CORRECTION FOR THE SATELLITES. COMPUTE THE FREQUENCIES

* Look for empty filenumber

      do unit=99,7,-1
	 inquire (unit=unit,opened=opened)
	 if (.not.opened) goto 10
      enddo
10    l=index(directory,' ')
      filenm=directory
      filenm(l:)='/tides.inp'
      open (unit,file=filenm,status='old',form='formatted')
      read (unit,1) s0,h0,p0,enp0,pp0,ds,dh,dp,dnp,dpp
      dtau = 365+dh-ds
      do 101 k = 1,ncomp
       read (unit,2) ii(k),jj(k),kk(k),ll(k),mm(k),nn(k),
     |            semi(k),nj(k)
       freq(k) = (ii(k)*dtau+jj(k)*ds+kk(k)*dh+ll(k)*dp+mm(k)*dnp+
     |                                      nn(k)*dpp)/365
       if (nj(k) .eq. 0) go to 101 
       read (unit,3) (ldel(k,j),mdel(k,j),ndel(k,j),ph(k,j),ee(k,j),
     |             ir(k,j),j=1,nj(k))
101   continue
      close (unit)
      return

555   continue   

* SET UP THE ARRAYS V(K),U(K),F(K) FOR THE COMPONENTS
      years=(mjd-mjd1976)/365
      s = s0+years*ds
      h = h0+years*dh
      p = p0+years*dp
      enp = enp0+years*dnp
      pp = pp0+years*dpp
      tau = mod(mjd,1d0)+h-s
      do k = 1,ncomp
         vdbl = ii(k)*tau+jj(k)*s+kk(k)*h+ll(k)*p+mm(k)*enp+
     |                                      nn(k)*pp+semi(k)
         iv = int(vdbl)
         iv = (iv/2)*2
         v(k) = vdbl-iv
         do j = 1,nj(k)
            uudbl = ldel(k,j)*p+mdel(k,j)*enp+ndel(k,j)*pp+ph(k,j)
            uu = uudbl-int(uudbl)
            uucos(k,j) = cos(twopi*uu)
            uusin(k,j) = sin(twopi*uu)
         enddo
      enddo

666   continue

*  SET ARRAYS U AND F (THESE ARE LATITUDE DEPENDENDENT)
      if (abs(xlat) .lt. 1d-5) then
         slat = xlat/abs(xlat)*1d-5
      else
         slat = xlat
      endif
      slat = sin(rad*slat)
      do k = 1,ncomp
         sumc = 1
         sums = 0
         do j = 1,nj(k)
            rr = ee(k,j)
            if (ir(k,j) .eq. 1) then
               rr = ee(k,j)*0.36309d0*(1-5*slat*slat)/slat
            else if (ir(k,j) .eq. 2) then
               rr = ee(k,j)*2.59808d0*slat
            endif
            sumc = sumc+rr*uucos(k,j)
            sums = sums+rr*uusin(k,j)
         enddo
         f(k) = sqrt(sumc*sumc+sums*sums)
         u(k) = atan2(sums,sumc)/twopi
      enddo
      return
      end

      function tidht2 (mjd,xlat,amp,phase)

*  COMPUTES THE P.O.L TIDAL CORRECTIONS AT THE DATA POINTS.
*  BASED ON A PROGRAM BY MIKE FOREMAN, INSTITUTE OF OCEAN SCIENCES,
*  CANADA

      include 'flather.inc'

      real*8 tidht2,mjd,mjdm,xlat,amp(ncomp),phase(ncomp)

      integer*4 icall
      integer*4 k,iml,imm
      save 

*  READ IN FROM FILE tides.inp AND SET FREQUENCIES ON FIRST CALL
      if (icall .eq. 0) then
         call astro2 (0,0d0,0d0)
         iml=-99999
         icall=1
      endif

*  CALCULATE V,U & F IF THE MONTH CHANGES, ELSE ONLY U & F
      imm=nint(mjd/30)
      if (imm.ne.iml) then
	 mjdm=imm*30
         iml = imm
         call astro2 (1,mjdm,0d0)
      endif
      call astro2 (2,mjdm,xlat)

*  CALCULATE THE HEIGHTS
      tidht2=0
      do k=1,ncomp
         tidht2=tidht2+f(k)*amp(k)*
     |        cos(twopi*(v(k)+(mjd-mjdm)*freq(k)+u(k)-phase(k)/360))
      enddo
      return
      end
