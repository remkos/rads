**FESTIDE -- Compute tides according to FES model
*+
      SUBROUTINE FESTIDE (UTC, LAT, LON, TIDE, TIDE_LP, TIDE_LOAD)
      REAL*8    UTC, LAT, LON, TIDE, TIDE_LP, TIDE_LOAD

* This routine makes the tidal predictions of ocean and load tide
* (optional) based on one of the FES models. This routine is heavily
* based on the routines by J.M. Molines and F. Lefevre for the FES99,
* FES2002 and FES2004 models.
*
* The input grids can be found in $ALTIM/data/$TIDE. The grids are
* in NetCDF format and were converted using the program "fes2nc".
*
* To initialize the computation, the function FESINIT should be
* called first. It allocates the appropriate amount of memory and
* loads the grids into memory. To release the memory for further
* use, call DALLOCF.
*
* Longitude and latitude are to be specified in degrees; time in UTC
* seconds since 1 Jan 1985. All predicted tides are output in meters.
* If the tide is requested in a point where it is not defined, NAN
* (Not-a-Number) is returned.
*
* Input arguments:
*  UTC      : UTC time in seconds since 1 Jan 1985
*  LAT      : Latitude (degrees)
*  LON      : Longitude (degrees)
*
* Output arguments:
*  TIDE     : Predicted short-period tide (m)
*  TIDE_LP  : Predicted long-period tide (m)
*  TIDE_LOAD: Predicted loading effect (m)
*-
* $Log: festide.f,v $
* Revision 1.23  2007/04/02 17:30:50  rads
* - Update to 2007 version of FES2004
*
* Revision 1.22  2006/09/29 15:08:50  rads
* - Send progress reports to standard error, not standard output
*
* Revision 1.21  2006/09/25 01:51:45  rads
* - Use NetCDF grids as input
* - Fix amplitude of M4 tide (was about 9% too small)
* - Use new LPETIDE (more waves)
*
* Revision 1.18  2006/08/07 17:24:59  rads
* - Removing obsolete Fortran code so it can be used with gfortran
*
* Revision 1.15  2005/10/27 20:14:43  rads
* - Added S1 and M4 short period tides and Mm, Mf, Mtm and MSqm long period tides
*
* Revision 1.13  2005/05/31 20:52:45  rads
* - Do all internal computation in double precision
*
* Revision 1.12  2005/03/10 10:57:44  eelco
* Removed calls to perror and added makefiles for compatibility with Intel Fortran Compiler 8.1
*
* Revision 1.11  2004/11/23 01:54:11  remko
* - Implementation supports FES2004
* - Automatically allocates memory, using FESINIT and FESFREE
* - Change of arguments in FESTIDE
*
* 10-Jun-2003 - Prepared for FES2002 (P1 component from grids);
*               switched M11/M12 admittance amplitudes;
*               fixed M12 frequency; allow grids with various scales.
*  9-Oct-2001 - Adapted by Remko Scharroo for DEOS
* 09-May-2001 - FES99 version by Fabien Lefevre (CLS)
* 16-Jun-1995 - First version by J. M. Molines (Grenoble)
*-----------------------------------------------------------------------
      include "festide.inc"
      include "nan.inc"
      real*8    mjd85,slon,delta,delta_max,t,lpetide
      parameter	(mjd85=46066d0,delta_max=24d0)
      complex*16 val(nw),w
      integer*4 kw,istat,mode

* When called with invalid time, reset time reference and return

      if (utc.gt.1d20) then
         t_nodal=1d30
	 return
      endif

* Compute long period tide. Use standard CTE waves for older versions
* of FES. When the FES model contains grids for long-period tides, use the
* all waves, except the monthly and fortnightly waves (Mm, Mf, Mtm and Msqm)

      mode=0
      if (nw_long.gt.0) mode=3
      tide_lp=lpetide(utc/86400d0+mjd85,lat,mode)

* If latitude out of range, bail out

      if (lat.lt.ymin .or. lat.gt.ymax) then
         tide=nan
	 tide_load=nan
	 return
      endif

* Limit longitude to interval of grids

      slon=lon
      if (slon.gt.xmax) then
	 slon=slon-360d0
      else if (slon.lt.xmin) then
	 slon=slon+360d0
      endif

* Conversion sec85 to days since 1900

      t=utc/86400d0+31046d0

* Determine time since last call to fes_astron (in hours).
* If time lapse is larger than delta_max, call astronimics again.

      delta=(t-t_nodal)*24d0
      if(abs(delta).ge.delta_max) then
        call fes_astron(t/36525d0)
	t_nodal=t
        delta=0d0
      endif

* Compute all ocean tide components

      tide=0d0

      call fes_interp(slon,lat,nw_prime+nw_extra+nw_long,
     |	val,tmp_s(ptr_ocean),dz(1,1),z0(1,1),nx,ny,dx,dy,istat)
      if (istat.eq.0) then
	tide=nan
	tide_lp=nan
      else
        do kw=1,nw
	  w = exp(cmplx(0d0,-(freq(kw)*delta+v0_u(kw))))
	  if (kw.ge.12 .and. kw.le.15) then
	     tide_lp = tide_lp + f(kw)*dble(w*val(kw))
	  else
             tide = tide + f(kw)*dble(w*val(kw))
	  endif
        enddo
      endif

* Compute all load tide components (if requested)

      tide_load=0d0
      if (.not.haveload) return

      call fes_interp(slon,lat,nw_prime,
     |	val,tmp_s(ptr_load),dz(1,2),z0(1,2),px,py,fx,fy,istat)
      if (istat.eq.0) then
	tide_load=nan
      else
        do kw=1,nw
	  w = exp(cmplx(0d0,-(freq(kw)*delta+v0_u(kw))))
          tide_load = tide_load + f(kw)*dble(w*val(kw))
        enddo
      endif
      end

*&FESINIT -- Initialize FES tide model
*+
      FUNCTION FESINIT (NAME, WANTLOAD)
      CHARACTER*(*) NAME
      LOGICAL	WANTLOAD
      INTEGER*4 FESINIT
*
* Allocate memory for FES tide modeling and read grids into memory.
* When WANTLOAD is .TRUE., loading tide grids are loaded and load tide
* will be computed. When .FALSE., FESTIDE will return a zero load tide.
* 
* Input arguments:
*  NAME     : Name of the FES tide model (FES95.2.1, FES99, or
*             FES2002, or FES2004)
*  WANTLOAD : Specify that load tide has to be computed
*
* Returned value:
*  FESINIT  : Pointer to memory buffer
*=
      include "festide.inc"
      integer	mallocf
      character*256 pathname
      integer*4	lex,lmod,lwa,lnblnk,kw,i,ptr,memloc
      logical	fes_read_tide
*
* These values are internal to the program and to the nodal correction
* routines. Do not change under any circumstances.
* Primary tides (nw_prime)
      data (wave(i),freq(i),i=1,9) /
     |     'Q1'     ,13.39866087990d0,
     |     'O1'     ,13.94303558000d0,
     |     'K1'     ,15.04106864000d0,
     |     '2N2'    ,27.89535481990d0,
     |     'N2'     ,28.43972952010d0,
     |     'M2'     ,28.98410422000d0,
     |     'K2'     ,30.08213728000d0,
     |     'S2'     ,30.00000000000d0,
     |     'P1'     ,14.95893136000d0/

* Additional small tides without loading components (nw_extra)
      data (wave(i),freq(i),i=10,11) /
     |     'M4'     ,57.9682084d0,
     |     'S1'     ,15.0000000d0/

* Long-period tides (nw_long)
      data (wave(i),freq(i),i=12,15) /
     |     'Mf'     ,1.09803310d0,
     |     'Mm'     ,0.54437470d0,
     |     'Mtm'    ,1.64240780d0,
     |     'MSqm'   ,2.11392880d0/

* Extra diurnal and semi diurnal are infered by admittance (nw_admit)
      data (wave(i),freq(i),i=16,33) /
     |     'Nu2'    ,28.51258314000d0,
     |     'Mu2'    ,27.96820844000d0,
     |     'L2'     ,29.52847892000d0,
     |     'T2'     ,29.95893332010d0,
     |     'Eps2'   ,27.4238337d0,
     |     'Lambda2',29.4556253d0,
     |     'Eta2'   ,30.6265120d0,
     |     '2Q1'    ,12.8542862d0,
     |     'Sigma1' ,12.9271398d0,
     |     'Ro1'    ,13.4715145d0,
     |     'M11'    ,14.4966939d0,
     |     'M12'    ,14.4874103d0,
     |     'Ki1'    ,14.5695476d0,
     |     'Pi1'    ,14.9178647d0,
     |     'Phi1'   ,15.1232059d0,
     |     'Teta1'  ,15.5125897d0,
     |     'J1'     ,15.5854433d0,
     |     'OO1'    ,16.1391017d0/

* Determine space needed to load tide grids
* nbytes = nwaves * nx * ny * 2 (im and re) * 2 (bytes)

      FESINIT=0
      nw_prime=9
      nw_extra=0
      nw_long=0
      nb_ocean=1441*721
      nb_load=nb_ocean
      if (name(:5).eq.'FES95') then
	 nw_prime=8
         nb_ocean=721*361
	 nb_load=nb_ocean
      else if (name(:5).eq.'FES99') then
	 nw_prime=8
      else if (name(:7).eq.'FES2002') then
      else if (name.eq.'FES2004') then
	 nw_extra=2
	 nw_long=4
         nb_ocean=2881*1441
	 nb_load=nb_ocean
      else
         nb_ocean=2881*1441
      endif
      nb_ocean=nb_ocean*(nw_prime+nw_extra+nw_long)*2*2
      nb_load=nb_load*nw_prime*2*2
      if (.not.wantload) nb_load=0

      if (mallocf(nb_ocean+nb_load,ptr_ocean).ne.0) then
         write (0,*) 'FESINIT: not able to allocate memory'
	 return
      endif
      FESINIT=ptr_ocean
      ptr_ocean=(ptr_ocean-memloc(tmp_s))/2

************************************************************************
* Read the data files for each tidal constituents. The path to the
* directory where the data files are located is given by $ALTIM/data,
* where $ALTIM is an environment valiable.
* All files are named model/{ocean,load}/wavename_{im,re}.grd
* where model is FES95.2.1, FES99 or FES2002 and wavename is M2,K1, etc.

      haveload=wantload
      lmod=lnblnk(name)
      pathname='/user/altim'
      call checkenv('ALTIM',pathname,lwa)
      pathname(lwa+1:)='/data/'//name(:lmod)//'/'
      lwa=lnblnk(pathname)

600   format ('(Loading FES tide: ',a,a,$)
601   format (a,$)
602   format (a,a,a)

* Eight or nine (nw_prime) primary waves are read from grids.
* Another 2 additional grids and 4 long period waves are also loaded.

      ptr=ptr_ocean
      write (0,600) pathname(:lwa),'{'
      do kw=1,nw_prime+nw_extra+nw_long
	lex=lnblnk(wave(kw))
	pathname(lwa+1:)=wave(kw)(:lex)//'_fes'//name(4:lmod)//'.nc'
	if (fes_read_tide(pathname,tmp_s(ptr),dz(kw,1),z0(kw,1),nx,ny,dx,dy))
     |		then
	   if (kw.gt.1) write (0,601) ','
	   write (0,601) wave(kw)(:lex)
	   ptr=ptr+2*nx*ny
	else
	   stop 'Error opening file. Fatal error.'
	endif
      enddo
      write (0,602) '}_fes',name(4:lmod),'.nc)'

* Read loading grids (when available)

      if (.not.haveload) goto 101
      write (0,600) pathname(:lwa),'{'
      do kw=1,nw_prime
        lex=lnblnk(wave(kw))
	pathname(lwa+1:)=wave(kw)(:lex)//'_drfes'//name(4:lmod)//'.nc'
	if (fes_read_tide(pathname,tmp_s(ptr),dz(kw,2),z0(kw,2),px,py,fx,fy))
     |		then
	   if (kw.gt.1) write (0,601) ','
	   write (0,601) wave(kw)(:lex)
	   ptr=ptr+2*px*py
        else
	   haveload=.false.
           goto 101
        endif
      enddo
      write (0,602) '}_drfes',name(4:lmod),'.nc)'

101   do kw=1,nw
        freq(kw)=freq(kw)*rad
      enddo

      t_nodal=1d30
      ptr_load=ptr_ocean+nb_ocean/2
      end

************************************************************************
* Perform bi-linear interpolation at point x,y from the gridded files.
* istat returns the number of points used for the interpolation

      subroutine fes_interp(x,y,ng,val,work,fact,offs,kx,ky,sx,sy,istat)
      include "festide.inc"
      complex*16 val(*)
      integer*4 istat,kw,i0,j0,i,j,kx,ky,ng
      integer*2 work(kx,ky,2,*)
      real*8	x,y,ptot,xij,yij,pds,fact(*),offs(*),sx,sy

      istat=0

      do kw=1,nw
        val(kw)=(0d0,0d0)
      enddo

      xij=(x-xmin)/sx+1
      yij=(y-ymin)/sy+1
      i0=min(kx-1,int(xij))
      j0=min(ky-1,int(yij))

      ptot=0d0

      do i=i0,i0+1
        do j=j0,j0+1
* Check if all major constituents are available:
          do kw=1,ng
            if(work(i,j,1,kw).eq.undef) goto 200
          enddo

          istat=istat+1
          pds=(1-abs(i-xij))*(1-abs(j-yij))
          ptot=ptot+pds

          do kw=1,ng
	    val(kw)=val(kw)+
     |		(work(i,j,1,kw)*fact(kw)+offs(kw))*exp(radj*work(i,j,2,kw))*pds
          enddo
200	  continue
        enddo
      enddo

      if (istat.eq.0 .or. ptot.eq.0d0) then
        istat=0
	return
      endif
      do kw=1,ng
        val(kw)=val(kw)/ptot
      enddo

* Infer additional constituents by admittance

* SEMI-DIURNAL (from Grenoble to take advantage of 2N2)

      val(16)= -0.0061047d0*val(7)+0.1568788d0*val(5)+0.0067557d0*val(6)	! Nu2 from N2,M2,K2
      val(17)=  0.0694400d0*val(7)+0.3515356d0*val(5)-0.0462783d0*val(6)	! Mu2 from N2,M2,K2
      val(18)=  0.0771378d0*val(7)-0.0516535d0*val(5)+0.0278699d0*val(6)	! L2 from N2,M2,K2
      val(19)=  0.1804802d0*val(7)-0.0201012d0*val(5)+0.0083315d0*val(6)	! T2 from N2,M2,K2
      val(20)=  0.53285d0  *val(4)-0.03304d0  *val(5) 				! Eps2 from 2N2,N2
      val(21)=  0.0165036d0*val(7)-0.0133078d0*val(5)+0.0077534d0*val(6)	! Lambda2 from N2,M2,K2
      val(22)= -0.0034925d0*val(6)+0.0831707d0*val(7)				! Eta2 from M2,K2

* DIURNALS (from Richard Ray perth2 program. Thank You, Richard!)

      val(23)=  0.263d0 *val(1)-0.0252d0*val(2)					! 2Q1 from Q1,O1
      val(24)=  0.297d0 *val(1)-0.0264d0*val(2)					! Sigma1 from Q1,O1
      val(25)=  0.164d0 *val(1)+0.0048d0*val(2)					! Ro1 from Q1,O1
      val(26)=  0.0389d0*val(2)+0.0282d0*val(3)					! M11 from O1,K1
      val(27)=  0.0140d0*val(2)+0.0101d0*val(3)					! M12 from O1,K1
      val(28)=  0.0064d0*val(2)+0.0060d0*val(3)					! K1 from O1,K1
      val(29)=  0.0030d0*val(2)+0.0171d0*val(3)					! Pi1 from O1,K1
      val(30)= -0.0015d0*val(2)+0.0152d0*val(3)					! Phi1 from O1,K1
      val(31)= -0.0065d0*val(2)+0.0155d0*val(3)					! Teta1 from O1,K1
      val(32)= -0.0389d0*val(2)+0.0836d0*val(3)					! J1 from O1,K1
      val(33)= -0.0431d0*val(2)+0.0613d0*val(3)					! OO1 from O1,K1

* P1 from Grenoble admittance code when grid is not provided

      if (nw_prime.lt.9)
     |val( 9)= -0.2387302d0*val(1)+0.1038608d0*val(2)+0.2892755d0*val(3)	! P1 from Q1,O1,K1

      end

************************************************************************
* Initialize some astronomic data and compute nodal corrections.

      subroutine fes_astron(tj)
      include "festide.inc"
      real*8	n,p,s,p1,pp,nu,xi,tt,nuprim,nusec,hp,r,iang,x1ra,hpi
      real*8	tj,tgn2,at1,at2,u,tgi2,coshn,sinhn,sin1n,sin2n
      integer*4 i

      tt=mod(180d0         +  360d0*36525d0*tj, 360d0)*rad	! Mean solar angle relative to Greenwich
      n =mod(259.1560563d0 - 1934.1423972d0*tj, 360d0)*rad	! Longitude of ascending lunar node
      hp=mod(280.1895015d0 +  36000.76892d0*tj, 360d0)*rad	! Mean solar longitude
      s =mod(277.0256206d0 +   481267.892d0*tj, 360d0)*rad	! Mean lunar longitude
      p1=mod(281.2208568d0 +     1.719175d0*tj, 360d0)*rad	! Longitude of solar perigee
      p =mod(334.3837214d0 + 4069.0322056d0*tj, 360d0)*rad	! Longitude of lunar perigee

      u=9.13694997d-1-3.5692561d-2*cos(n)
      iang=acos(u)
      tgn2=tan(n/2)
      at1=atan(1.01883d0*tgn2)
      at2=atan(6.4412d-1*tgn2)
      xi=-at1-at2+n
      if(n.gt.pi) xi=xi-2*pi
      nu=at1-at2

* For constituents L2, K1, K2

      tgi2=tan(iang/2)
      pp=p-xi
      x1ra=sqrt(1-12*tgi2**2*cos(2*pp)+36*tgi2**4)
      r=atan(sin(2*pp)/(1/(6*tgi2**2)-cos(2*pp)))
      nuprim=atan(sin(2*iang)*sin(nu)/
     |     (sin(2*iang)*cos(nu)+3.347d-1))
      nusec=0.5d0*atan(((sin(iang)**2)*sin(2*nu))
     |     /(sin(iang)**2*cos(2*nu)+7.27d-2))

* Compute nodal corrections from Schureman (1958)

      sinhn=sin(iang/2d0)
      sin1n=sin(iang)
      sin2n=sin(2d0*iang)
      coshn=cos(iang/2d0)
      f( 1)=sin1n*coshn**2/0.38d0
      f( 2)=f( 1)
      f( 3)=sqrt(0.8965d0*sin2n**2+0.6001d0*sin2n*cos(nu)+0.1006d0)
      f( 4)=coshn**4/0.9154d0
      f( 5)=f( 4)
      f( 6)=f( 4)
      f( 7)=sqrt(19.0444d0*sin1n**4+2.7702d0*sin1n**2*cos(2*nu)+.0981d0)
      f( 8)=1d0
      f( 9)=f( 8)
      f(10)=f( 4)**2
      f(11)=1d0
      f(12)=sin1n**2/0.1578d0
      f(13)=(2d0/3d0-sin1n**2)/0.5021d0
      f(14)=f(12)
      f(15)=f(12)
      f(16)=f( 4)
      f(17)=f( 4)
      f(18)=f( 4)*x1ra
      f(19)=f( 8)
      f(20)=f( 8)
      f(21)=f( 4)
      f(22)=sin1n**2/0.1565d0
      f(23)=f( 1)
      f(24)=f( 1)
      f(25)=f( 1)
      f(26)=sin2n/0.7214d0
      f(27)=f( 1)
      f(28)=f(26)
      f(29)=f( 8)
      f(30)=f( 8)
      f(31)=f(26)
      f(32)=f(26)
      f(33)=sin1n*sinhn*sinhn/0.0164d0

* Compute V0+u from Schureman (1958)

      hpi=pi/2
      v0_u( 1)=tt-3*s+hp+p+hpi+2*xi-nu		! Q1
      v0_u( 2)=tt-2*s+hp+hpi+2*xi-nu		! O1
      v0_u( 3)=tt+hp-hpi-nuprim			! K1
      v0_u( 4)=2*(tt-2*s+hp+p+xi-nu)		! 2N2
      v0_u( 5)=2*tt-3*s+2*hp+p+2*xi-2*nu	! N2
      v0_u( 6)=2*(tt-s+hp+xi-nu)		! M2
      v0_u( 7)=2*(tt+hp-nusec)			! K2
      v0_u( 8)=2*tt				! S2
      v0_u( 9)=tt-hp+hpi			! P1
      v0_u(10)=4*(tt-s+hp+xi-nu)		! M4
      v0_u(11)=tt				! S1
      v0_u(12)=2*(s-xi)				! Mf
      v0_u(13)=s-p				! Mm
      v0_u(14)=3*s-p-2*xi			! Mtm
      v0_u(15)=2*(2*s-hp-xi)			! MSqm
      v0_u(16)=2*tt-3*s+4*hp-p+2*xi-2*nu	! Nu2
      v0_u(17)=2*tt-4*s+4*hp+2*xi-2*nu		! Mu2
      v0_u(18)=2*tt-s+2*hp-p+pi+2*xi-2*nu-r	! L2
      v0_u(19)=2*tt-hp+p1			! T2
      v0_u(20)=2*tt-5*s+4*hp+p			! Eps2
      v0_u(21)=2*tt-s+p+pi+2*xi-2*nu		! Lambda2
      v0_u(22)=2*tt+s+2*hp-p-2*nu		! Eta2
      v0_u(23)=tt-4*s+hp+2*p+hpi+2*xi-nu	! 2Q1
      v0_u(24)=tt-4*s+3*hp+hpi+2*xi-nu		! Sigma1
      v0_u(25)=tt-3*s+3*hp-p+hpi+2*xi-nu	! Ro1
      v0_u(26)=tt-s+hp+p-hpi-nu			! M11
      v0_u(27)=tt-s+hp-p-hpi+2*xi-nu		! M12
      v0_u(28)=tt-s+3*hp-p-hpi-nu		! Ki1
      v0_u(29)=tt-2*hp+p1+hpi			! Pi1
      v0_u(30)=tt+3*hp-hpi			! Phi1
      v0_u(31)=tt+s-hp+p-hpi-nu			! Teta1
      v0_u(32)=tt+s+hp-p-hpi-nu			! J1
      v0_u(33)=tt+2*s+hp-hpi-2*xi-nu		! OO1

      do i=1,nw
         v0_u(i)=mod(v0_u(i),2*pi)
      enddo
      end
************************************************************************
      function fes_read_tide(pathname,work,fact,offs,kx,ky,sx,sy)

* Reads real and imaginary grids in NetCDF format

      include "festide.inc"
      include "netcdf.inc"
      integer*2 work(*)
      integer*4 kx,ky,ncid
      character*(*) pathname
      logical	fes_read_tide
      real*8	fact,offs,sx,sy,dummy(2)

      fes_read_tide=.false.

* Open and read NetCDF grid of amplitude and phase

      if (nf_open(pathname,nf_nowrite,ncid).ne.nf_noerr) return
      if (nf_inq_dimlen(ncid,1,kx).ne.nf_noerr) return
      if (nf_inq_dimlen(ncid,2,ky).ne.nf_noerr) return
      if (nf_get_att_double(ncid,1,"actual_range",dummy).ne.nf_noerr)
     |	return
      xmin=dummy(1) ; xmax=dummy(2)
      if (nf_get_att_double(ncid,2,"actual_range",dummy).ne.nf_noerr)
     |	return
      ymin=dummy(1) ; ymax=dummy(2)
      if (nf_get_att_double(ncid,3,"scale_factor",fact).ne.nf_noerr)
     |	fact=1d0
      if (nf_get_att_double(ncid,3,"add_offset",offs).ne.nf_noerr)
     |	offs=0d0
      if (nf_get_var_int2(ncid,3,work(1)).ne.nf_noerr) return
      if (nf_get_var_int2(ncid,4,work(kx*ky+1)).ne.nf_noerr) return
      if (nf_close(ncid).ne.nf_noerr) return

      sx=(xmax-xmin)/(kx-1)
      sy=(ymax-ymin)/(ky-1)
      fes_read_tide=.true.
      end
