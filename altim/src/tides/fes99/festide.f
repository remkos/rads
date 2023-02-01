c######################################################################
c  Main routines for the FES prediction software.
c
c  Version 3.2  : May, 9 2001
c                 update real*4, real*8, integer*4 to be compatible
c                 with most of Unix and Linux platform 
c  Minor changes: revision by Fabien LEFEVRE
c 
c  Version 3.1  : November, 1 2000
c  Major changes: full revision by Fabien LEFEVRE
c 
c  Files are distributed as ASCII files, but can be converted to binary
c  files by the program ascii2bin. Doing so reduces both the required disk
c  space and system time for reading the files.
c  Depending on your choice, modify the iascii variable
c  in the subroutine fes_tide
c 
c  This software have been tested on SUN Solaris 2.7
c  It is provided without any guarantees.
c 
c  For bug reports, please contact :
c  ---------------------------------
c  Fabien LEFEVRE 
c
c  CLS
c  http://www.cls.fr
c  Direction Océanographie Spatiale
c  8-10, rue Hermès - Parc Technologique du Canal
c  31526 Ramonville Saint-Agne cedex - France
c  Tel: +33 (0)5 61 39 37 45 Fax: +33 (0)5 61 39 37 82
c  e-mail: Fabien.Lefevre@cls.fr
c
c  NOTE: This software is based on the former versions
c        developed by Jean-Marc MOLINES and Florent LYARD
c ######################################################################

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c    SUBROUTINES
c    
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        


c########################################################################
c      
      subroutine fes_tide(model_name,wave_path_name,
     &     tide,tide_lp,tideload,lat,lon,jnasa,heure,iascii,istat)
c
c------------------------------------------------------------------------
c  ROUTINE : fes_tide
c
c  DESCRIPTION : tide prediction routine (end user routine):
c          tide      = predicted tide in cm (output)
c          tide_lp   = predicted long period tide in cm (output)
c          tide_load = predicted loading effect in cm (output)
c          lat       = latitude (between -90 et +90 deg ) (input)
c          lon       = longitude (either between -180 180 or 0 360).
c          jnasa     = Day for prediction (as in TP GDR )
c          heure     = Hour during the day jnasa
c          istat     = Number of valid points in the interpolation
c
c  PROGRAMMERS : J.M. MOLINES and F. LEFEVRE
c 
c  DATE : 16/06/95 - 09/05/01
c------------------------------------------------------------------------
c Conversion date :
c - JD  = Julian Day
c - MJD = Modified Julian Day : MJD = JD - 2400000.5
c
c NASA day reference : 1 jan 1958 (JD = 2436204.5, MJD = 36204)
c CNES day reference : 1 jan 1985 (JD = 2446066.5, MJD = 46066)
c------------------------------------------------------------------------
      include 'common.h' 
c
      real*8        t1,t2,heure
      character*(*) model_name,wave_path_name
      real*4        lat,lon,tide,tide_lp,tideload
      real*8        dlat, dtide_lp
      integer*4     jnasa, istat, ifirst, iascii

      data ifirst /0/
      
      iascbin=iascii

c-- Read the data files at first call. This is quite long
      if (ifirst.eq.0) then
        model=model_name
        wave_path=wave_path_name
        call init_data
        ifirst=1
      endif

      if (lon.ge.180.) lon=lon-360.

c-- Conversion TOPEX date to adequate date for prediction
      call convert_date(jnasa,heure,t1)
c      print*, 'convert_date ',jnasa,heure,t1

c-- Compute diurnal and semi-diurnal tide
      call tidal_comp(t1,lon,lat,tide,tideload,istat)

c-- Compute long period tide
      t2=(((jnasa+36204)*1.d0)*24.d0+heure)*3600.d0
      dlat=dble(lat)
      call lpeqmt(t2, dlat, dtide_lp)
      tide_lp=dtide_lp

      return
      end


c########################################################################
c      
      subroutine astronomics(tj)
c
c------------------------------------------------------------------------
c  ROUTINE : astronomics
c
c  DESCRIPTION : This program initialize some astronomic data useful for
c                nodal corrections.
c
c  PROGRAMMERS : J.M. MOLINES and F. LYARD
c 
c  DATE : 19/01/94
c------------------------------------------------------------------------
      include 'common.h'          
      
      real*8 tj,rad,tgn2,at1,at2,u,tgi2
      real*8 ct0,ct1
      real*8 cn0,cn1
      real*8 cs0,cs1
      real*8 ch0,ch1
      real*8 cps0,cps1
      real*8 cp0,cp1  

c-- tt mean solar angle relative to Greenwich -------------------------
      
      ct0=180.D+00
      ct1=360.D+00*3.6525D+04
      tt=dmod(ct0+ct1*tj,360.D+00)

c-- hp longitude of ascending lunar node ------------------------------
      
      cn0=  259.1560563D+00
      cn1= -1934.1423972D+00   
      n=dmod(cn0+cn1*tj,360.D+00)

c-- hp mean solar longitude -------------------------------------------

      ch0= 280.1895015D+00
      ch1= 36000.76892D+00 
      hp=dmod(ch0+ch1*tj,360.D+00)
      
c-- s mean lunar longitude --------------------------------------------

      cs0= 277.0256206D+00
      cs1= 481267.892D+00
      s=dmod(cs0+cs1*tj,360.D+00)

c-- p1 longitude of solar perigee -------------------------------------
      
      cps0=281.2208568D+00
      cps1=1.719175D+00   
      p1=dmod(cps0+cps1*tj,360.D+00)

c-- p longitude of lunar perigee --------------------------------------

      cp0=334.3837214D+00
      cp1=4069.0322056D+00
      p=dmod(cp0+cp1*tj,360.D+00)

      n=rad(n)
      hp=rad(hp)
      s=rad(s)
      p1=rad(p1)
      p=rad(p)

      u=9.13694997D-01-3.5692561D-02*dcos(n)
      iang=dacos(u)
      tgn2=dtan(n/2.)
      at1=datan(1.01883D+00*tgn2)
      at2=datan(6.4412D-01*tgn2)
      xi=-at1-at2+n
      if(n.gt.pi) then
        xi=xi-2.D+00*pi
      endif
      nu=at1-at2
      
c-- for constituents l2,k1,k2 -----------------------------------------

      tgi2=dtan(iang/2.D+00)     
      pp=p-xi
      x1ra=dsqrt(1.-12.D+00*tgi2**2*cos(2.D+00*pp)+36.D+00*tgi2**4)
      r=datan(dsin(2.*pp)/(1./(6.*tgi2**2)-dcos(2.*pp)))
      nuprim=datan(dsin(2.D+00*iang)*dsin(nu)/
     &     (dsin(2.D+00*iang)*dcos(nu)+3.347D-01)) 
      nusec=0.5*datan(((dsin(iang)**2.D+00)*dsin(2.D+00*nu))
     &     /(dsin(iang)**2.D+00*dcos(2.*nu)+7.27D-02))                 

      return
      end


c########################################################################
c      
      subroutine calend1(njd,nd,nm,na)
c
c------------------------------------------------------------------------
c  ROUTINE : calend1
c
c  DESCRIPTION : Conversion from Jcnes to dd/mm/aaaa
c
c  PROGRAMMERS : Origine GRGS Toulouse, C. BROSSIER
c 
c  DATE : 19/01/94
c------------------------------------------------------------------------
      integer*4 njd, nd, nm, na, n, njul, nj, nb, j, nm1, m, ndj, nj3
      dimension n(12)
      data n /31,28,31,30,31,30,31,31,30,31,30,31/
      njul=njd+1
      na=njul/365
      nj=njul-na*365
      nb=(na+1)/ 4
      nj=nj-nb
      if(nj.gt.0) go to 1
      na=na+1949
      nm=12
      nd=nj+31
      return
    1 j=na-2-nb*4
      na=na+1950
      if(j.lt.0) go to 5000
 4000 if(60-nj)4500,7000,5000
 4500 nm1=60
      m=3
      go to 6000
 5000 nm1=0
      m=1
 6000 ndj=nm1+n(m)
      nj3=nj-ndj
      if(nj3.le.0) go to 8000
 6500 m=m+1
      nm1=ndj
      go to 6000
 7000 nm=2
      nd=29
      return
 8000 nm=m
      nd=nj-nm1
 9000 return
      end


c########################################################################
c      
      subroutine convert_date(jnasa,heure,t1)
c
c------------------------------------------------------------------------
c  ROUTINE : convert_date
c
c  DESCRIPTION : Convert from Jnasa to Julian century as used in
C                Schureman
c
c  PROGRAMMERS : J.M. MOLINES 
c 
c  DATE : 13/03/1995
c------------------------------------------------------------------------
c      
      integer*4 jnasa, ida, imonth, ian, iday, j
      integer*4 tojul
      real*8 heure,t1,julian
      real*4   sec
c
c CAUTION 2922 = Number of days between NASA and CNES ...
c
      call calend1(jnasa+2922,ida,imonth,ian)
      iday=tojul(ian,imonth,ida)
      j=iday-1000*(ian-1900)
      sec=heure*3600.d0
      t1=julian(ian,j,sec)
      return
      end


c########################################################################
c      
      subroutine convert_ymd2nasa(jnasa,iyear,imonth,iday)
c
c------------------------------------------------------------------------
c  ROUTINE : convert_ydm2nasa
c
c  DESCRIPTION : Convert from year/month/day to NASA day
c
c  PROGRAMMERS : F. LEFEVRE
c 
c  DATE : 06/11/00
c------------------------------------------------------------------------
c      
      integer*4 jnasa, imonth, iyear, iday, j, ida
      real*4   sec
      real*8 heure, t1
      integer*4 tojul
      real*8 julian

c   1 jan 1950 = 18262.0 julian days since 1 jan 1900
c   2922 = number of days between NASA et CNES ...
        ida=tojul(iyear,imonth,iday)
        j=ida-1000*(iyear-1900)
        heure=0.0
        sec=heure*3600.d0
        t1=julian(iyear,j,sec)

        jnasa=nint(t1*365.25d2)-18262-2922

      return
      end


c########################################################################
c      
      subroutine tidal_comp(t,rlon,rlat,vtide,vtideload,istat)
c
c------------------------------------------------------------------------
c  ROUTINE : tidal_comp
c
c  DESCRIPTION : driver for tidal computation.
c                The result is returned in tide.
c                tide is an array 
c
c  PROGRAMMERS : J.M. MOLINES and F. LYARD
c 
c  DATE : 13/13/95
c------------------------------------------------------------------------
      include 'common.h' 
      real*4 tide, tideload
      dimension tide(extended_max),  tideload(extended_max)        

c-----------------------------------------------------------------------
c istat = Number of grid points for data interpolation
c istat = 0  point (rlon,rlat) is out of the gridded area
c istat = 1,2,3  point (rlon,rlat) is on a frontier of the area
c istat = 4  point (rlon,rlat) is fully in the domain
c-----------------------------------------------------------------------

      real*4    rlon, rlat, vtide, vtideload
      integer*4 istat
      real*8  delta,t
      real*4    phi
      real*4    h, hload
      integer*4 nb   
      
      delta=(t-t_nodal)*3.6525d+04*24.D+00

      if(abs(delta).gt.delta_max) then
        call init_corrections(t)
        delta=0.D+00
      endif

c-----------------------------------------------------------------------
c      sreal(nb) = real*4 part of the nb th pure tide wave
c      simag(nb) = imaginary  part of the nb th pure tide wave
c      srealload(nb) = real*4 part of the nb th radial tide wave
c      simagload(nb) = imaginary  part of the nb th radial tide wave
c-----------------------------------------------------------------------

      call interpolation(rlon,rlat,sreal,simag,wrp,wip,istat)
      call interpolation(rlon,rlat,srealload,simagload,wrpload,
     & wipload,istat)    
      
      if (istat.ne.0) then
        h=0.0
        hload=0.0
        do nb=1,nbwave
          phi=real(freq(nb)*delta)+v0_u(nb)  
          phi=mod(phi,2.*pi)
          if (phi.lt.0.) phi=phi+2.*pi
          h=h+f(nb)*(sreal(nb)*cos(phi)+simag(nb)*sin(phi))
          tide(nb)=h
          hload=hload+f(nb)*(srealload(nb)*cos(phi)
     &          +simagload(nb)*sin(phi))
          tideload(nb)=hload
        enddo 
      else
        do nb=1,nbwave
          tide(nb)=spec
          tideload(nb)=spec
        enddo
      endif
      vtide=tide(nbwave)
      vtideload=tideload(nbwave)
    
      return
      end


c########################################################################
c      
      subroutine init_corrections(t0)
c
c------------------------------------------------------------------------
c  ROUTINE : init_corrections
c
c  DESCRIPTION : compute nodal corrections
c
c  PROGRAMMERS : J.M. MOLINES and F. LYARD
c 
c  DATE : 19/01/94
c------------------------------------------------------------------------
      include 'common.h'          
      
      real*8 t0
      integer*4 nb
      
      call nodalc(t0)     
      
      do nb=1,nbwave   
c-->    convert v0_u from degrees to radians
        v0_u(nb)=v0_u(nb)*pi/180.
      enddo   

      t_nodal=t0

      return
      end


c########################################################################
c      
      subroutine interpolation(xx,yy,sre,sim,wr,wi,istat)
c
c------------------------------------------------------------------------
c  ROUTINE : interpolation
c
c  DESCRIPTION : Perform bi-linear interpolation at point xx,yy from the
c                gridded files.
c                istat returns the number of points used for the
c                interpolation
c
c  August 2000 : no more constituents added by admittance 
c
c  PROGRAMMERS : J.M. MOLINES 
c 
c  DATE : 15/06/1995
c------------------------------------------------------------------------
      include 'common.h'   
      
      real*4 sre(extended_max), sim(extended_max)
      real*4 wr(1441,721,8), wi(1441,721,8)

c-- Variables for admittance computation
      real*4 aap1,bbp1,ccp1
      real*4 aamu2,aanu2,aal2,aat2,aalda2,bbmu2,bbnu2,bbl2,bbt2,bblda2
      real*4 ccmu2,ccnu2,ccl2,cct2,cclda2
      common /admitancesaa/aamu2,aanu2,aal2,aat2,aalda2,aap1
      common /admitancesbb/bbmu2,bbnu2,bbl2,bbt2,bblda2,bbp1
      common /admitancescc/ccmu2,ccnu2,ccl2,cct2,cclda2,ccp1
c-- End of Variables for admittance computation

      real*4    xx, yy, x, y, xmax, ptot, xij, yij, pds
      integer*4 istat, nb, i0, j0, i, j
      
      istat=0     

      if (nbwave.eq.0) return

      do nb=1,nbwave 
        sre(nb)=0.0
        sim(nb)=0.0  
      enddo           
      
      x=xx      
      y=yy   

      xmax=xmin+(ni-1)*dx
      if (x.lt.xmin) x=x+360.0
      if (x.gt.xmax) x=x-360.0

      i0=int((x-xmin)/dx)+1
      j0=int((y-ymin)/dy)+1 
      
      ptot=0.0

      do i=i0,i0+1 
        if(i.lt.0.or.i.gt.ni) goto 100
        do j=j0,j0+1 
          if(j.lt.1.or.j.gt.nj) goto 200
c Check if the 8 Major constituents are available:
         do nb=1,8
            if(wr(i,j,nb).gt.(spec-0.1).
     &           or.wi(i,j,nb).gt.(spec-0.1)) goto 200
          end do
c all minor constituents will be deduced from the major ones.
c
          istat=istat+1
          xij=xmin+(i-1)*dx  
          yij=ymin+(j-1)*dy
c this line modified for taking into account interpolation near the
C origin
          pds=(dx-mod(abs(x-xij),dx))*(dy-abs(y-yij))
          ptot=ptot+pds   

         do nb=1,8
            sre(nb)=sre(nb)+wr(i,j,nb)*pds    
            sim(nb)=sim(nb)+wi(i,j,nb)*pds
          enddo
 200      continue
        enddo 
 100    continue 
      enddo

      do nb=1,nbwave
        if(ptot.ne.0.0) then
          sre(nb)=sre(nb)/ptot  
          sim(nb)=sim(nb)/ptot  
        else    
          istat=0
          sre(nb)=spec  
          sim(nb)=spec  
        endif
      enddo 
c
      
c-- Variables for admittance computation
      if (istat.ne.0) then
c
c infer additional constituents by admittance
c
c DIURNALS (from Richard Ray perth2 program. Thank You, Richard!)
c---------
c  from Q1 and O1 (1-2)
c    2Q1
        sre(17) = 0.263 *sre(1) - 0.0252*sre(2)
        sim(17) = 0.263 *sim(1) - 0.0252*sim(2)
c    sigma1
        sre(18) = 0.297 *sre(1) - 0.0264*sre(2)
        sim(18) = 0.297 *sim(1) - 0.0264*sim(2)
c    rho1
        sre(19) = 0.164 *sre(1) + 0.0048*sre(2)
        sim(19) = 0.164 *sim(1) + 0.0048*sim(2)
c  from O1 and K1  (2-3)
c    M11
        sre(20) = 0.0140*sre(2) + 0.0101*sre(3)
        sim(20) = 0.0140*sim(2) + 0.0101*sim(3)
c    M12
        sre(21) = 0.0389*sre(2) + 0.0282*sre(3)
        sim(21) = 0.0389*sim(2) + 0.0282*sim(3)
c    chi1
        sre(22) = 0.0064*sre(2) + 0.0060*sre(3)
        sim(22) = 0.0064*sim(2) + 0.0060*sim(3)
c    pi1
        sre(23) = 0.0030*sre(2) + 0.0171*sre(3)
        sim(23) = 0.0030*sim(2) + 0.0171*sim(3)
c    phi1
        sre(24) =-0.0015*sre(2) + 0.0152*sre(3)
        sim(24) =-0.0015*sim(2) + 0.0152*sim(3)
c    theta1
        sre(25) =-0.0065*sre(2) + 0.0155*sre(3)
        sim(25) =-0.0065*sim(2) + 0.0155*sim(3)
c   J1
        sre(26) =-0.0389*sre(2) + 0.0836*sre(3)
        sim(26) =-0.0389*sim(2) + 0.0836*sim(3)
c   OO1
        sre(27) =-0.0431*sre(2) + 0.0613*sre(3)
        sim(27) =-0.0431*sim(2) + 0.0613*sim(3)
c
c    P1  from Grenoble admittance code
        sre(9) =  aap1*sre(1)+bbp1*sre(2)+ccp1*sre(3)
        sim(9) =  aap1*sim(1)+bbp1*sim(2)+ccp1*sim(3)
        
c SEMI-DIURNAL (from Grenoble to take advantage of 2N2)
c ------------
c  from 2N2 -N2 (4-5)
c    eps2 
        sre(14) = 0.53285 *sre(4) - 0.03304*sre(5)
        sim(14) = 0.53285 *sim(4) - 0.03304*sim(5)
c  from M2 - K2 (6-7)
c   eta2
        sre(16) = -0.0034925 *sre(6) + 0.0831707*sre(7)
        sim(16) = -0.0034925 *sim(6) + 0.0831707*sim(7)
c
c  from N2 -M2- K2 by spline admittances (see GRL 18(5):845-848,1991)
c
c   mu2
        sre(11) = aamu2*sre(7)+bbmu2*sre(5)+ccmu2*sre(6)
        sim(11) = aamu2*sim(7)+bbmu2*sim(5)+ccmu2*sim(6)
c   nu2
        sre(10) = aanu2*sre(7)+bbnu2*sre(5)+ccnu2*sre(6)
        sim(10) = aanu2*sim(7)+bbnu2*sim(5)+ccnu2*sim(6)
c   lda2
        sre(15) = aalda2*sre(7)+bblda2*sre(5)+cclda2*sre(6)
        sim(15) = aalda2*sim(7)+bblda2*sim(5)+cclda2*sim(6)
c   L2
        sre(12) = aal2*sre(7)+bbl2*sre(5)+ccl2*sre(6)
        sim(12) = aal2*sim(7)+bbl2*sim(5)+ccl2*sim(6)
c   T2
        sre(13) = aat2*sre(7)+bbt2*sre(5)+cct2*sre(6)
        sim(13) = aat2*sim(7)+bbt2*sim(5)+cct2*sim(6)

      endif
c-- End of Variables for admittance computation
      return
      end


c########################################################################
c      
      subroutine nodal_a      
c
c------------------------------------------------------------------------
c  ROUTINE : nodal_a
c
c  DESCRIPTION :Compute nodal corrections from SCHUREMAN (1958)
c     
c  CAUTION     : indexes used in this routine are internal to the code
c                and corresponds to the !original! ondes.dat file.
c      june'95 : new secondary constituents added.
c
c  PROGRAMMERS : J.M. MOLINES
c 
c  DATE : 15/06/1995
c------------------------------------------------------------------------
      include 'common.h'     
      
      integer*4 i

      do 100 i=1,nbwave 
        if(num(i).eq.1) goto 11 
        if(num(i).eq.2) goto 12 
        if(num(i).eq.3) goto 13  
        if(num(i).eq.5) goto 14
        if(num(i).eq.6) goto 14
        if(num(i).eq.7) goto 14
        if(num(i).eq.8) goto 14
        if(num(i).eq.9) goto 14
        if(num(i).eq.11) goto 15 
        if(num(i).eq.12) goto 12
        if(num(i).eq.13) goto 12 
        if(num(i).eq.14) goto 16 
        if(num(i).eq.27) goto 11
c
        if(num(i).eq.60) goto 12
        if(num(i).eq.61) goto 14
        if(num(i).eq.62) goto 14
        if(num(i).eq.63) goto 17
        if(num(i).eq.64) goto 17
        if(num(i).eq.65) goto 11
        if(num(i).eq.66) goto 11
        if(num(i).eq.67) goto 11
        if(num(i).eq.68) goto 18
        if(num(i).eq.69) goto 11
        if(num(i).eq.70) goto 18
        if(num(i).eq.71) goto 12
        if(num(i).eq.72) goto 12
        if(num(i).eq.73) goto 18
        if(num(i).eq.74) goto 18
        if(num(i).eq.75) goto 19
 11     f(i)=dsin(iang)*dcos(iang/2)**2/0.38  
        goto 120
 12     f(i)=1 
        goto 120 
 13     f(i)=sqrt(0.8965*dsin(2*iang)**2+0.6001*dsin(2*iang)
     &       *dcos(nu)+0.1006)          
        goto 120
 14     f(i)=real(dcos(iang/2.D+00)**4/9.154D-01)
        goto 120
 15     f(i)=dcos(iang/2)**4/0.9154*x1ra        
        goto 120 
 16     f(i)=sqrt(19.0444*dsin(iang)**4
     &       +2.7702*dsin(iang)**2*dcos(2*nu)+0.0981)
        goto 120
 17     f(i)=real(dsin(iang)*dsin(iang)/0.1565d0) 
        goto 120      
 18     f(i)=real(dsin(2*iang)/0.7214d0) 
        goto 120
 19     f(i)=real(dsin(iang)*dsin(iang/2.d0)*dsin(iang/2.d0)/0.0164d0) 
        goto 120
        
 120    continue   
 100  continue
      return
      end


c########################################################################
c      
      subroutine  nodal_G  
c
c------------------------------------------------------------------------
c  ROUTINE : nodal_g
c
c  DESCRIPTION : compute V0+u from Schureman (1958)
c       june 1995: new secondary constituents added
c
c
c  PROGRAMMERS : J.M. MOLINES and F. LYARD
c 
c  DATE : 16/06/95
c------------------------------------------------------------------------
      include 'common.h'          

      real*8 deg
      
      integer*4 i
      real*4 v0, u

      n=deg(n)
      hp=deg(hp)
      s=deg(s)    
      p1=deg(p1)
      p=deg(p)
      xi=deg(xi)
      nu=deg(nu)
      nuprim=deg(nuprim)
      nusec=deg(nusec)
      r=deg(r) 

      do 140 i=1,nbwave
c-------------------O1--------------------------------------------------
        if(num(i).eq.1) then
          v0=tt-2*s+hp+90.
          u=2*xi-nu
          v0_u(i)=tt-2*s+hp+90.+2*xi-nu
          goto 130
        endif

c-------------------P1--------------------------------------------------
        if(num(i).eq.2) then  
          v0=tt-hp+90.
          u=0
          v0_u(i)=tt-hp+90.
          goto 130
        endif      

c-------------------K1--------------------------------------------------
        if(num(i).eq.3) then           
          v0=tt+hp-90.                                             
          u=-nuprim
          v0_u(i)=tt+hp-90.-nuprim
          goto 130 
        endif
        
c-------------------2N2------------------------------------------------
        if(num(i).eq.5) then
          v0=2*tt-4*s+2*hp+2*p
          u=2*xi-2*nu
          v0_u(i)=2*tt-4*s+2*hp+2*p+2*xi-2*nu
          goto 130
        endif                                                           

c-------------------Mu2-------------------------------------------------
        if(num(i).eq.6) then
          v0=2*tt-4*s+4*hp
          u=2*xi-2*nu
          v0_u(i)=2*tt-4*s+4*hp+2*xi-2*nu
          goto 130 
        endif

c-------------------N2--------------------------------------------------
        if(num(i).eq.7) then
          v0=2*tt-3*s+2*hp+p
          u=2*xi-2*nu
          v0_u(i)=2*tt-3*s+2*hp+p+2*xi-2*nu
          goto 130
        endif

c-------------------Nu2-------------------------------------------------
        if(num(i).eq.8) then
          v0=2*tt-3*s+4*hp-p
          u=2*xi-2*nu
          v0_u(i)=2*tt-3*s+4*hp-p+2*xi-2*nu
          goto 130
        endif     

c-------------------M2--------------------------------------------------
        if(num(i).eq.9) then
          v0=2.*tt-2*s+2*hp
          u=2*xi-2*nu
          v0_u(i)=2.*tt+real(-2.D+00*s+2.D+00*hp+2.D+00*xi-2.D+00*nu)
          goto 130
        endif

c------------------L2---------------------------------------------------
        if(num(i).eq.11) then
          v0=2*tt-s+2*hp-p+180
          u=2*xi-2*nu-r
          v0_u(i)=2*tt-s+2*hp-p+180+2*xi-2*nu-r
          goto 130
        endif

c------------------T2---------------------------------------------------
        if(num(i).eq.12) then
          v0=2*tt-hp+p1
          u=0
          v0_u(i)=2*tt-hp+p1
          goto 130
        endif

c------------------s2---------------------------------------------------
        if(num(i).eq.13) then
          v0=2*tt
          u=0
          v0_u(i)=2*tt
          goto 130
        endif

c------------------K2---------------------------------------------------
        if(num(i).eq.14) then
          v0=2*tt+2*hp
          u=-2*nusec
          v0_u(i)=2*tt+2*hp-2*nusec
          goto 130
        endif

c------------------Q1---------------------------------------------------
        if(num(i).eq.27) then                 
          v0=tt-3*s+hp+p+90.
          u=2*xi-nu
          v0_u(i)=tt-3*s+hp+p+90.+2*xi-nu
          goto 130
        endif

c------------------Eps2-------------------------------------------------
        if(num(i).eq.60) then                 
          v0=2*tt-5*s+4*hp+p
          u=0.
          v0_u(i)=v0+u
          goto 130
        endif

c------------------Lambda2-------------------------------------------------
        if(num(i).eq.61) then                 
          v0=2*tt-s+p+180.
          u=2*xi-2*nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------l21--------------------------------------------------
        if(num(i).eq.62) then                 
          v0=2*tt-s +2*hp-p+180
          u=-2*nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------l22--------------------------------------------------
        if(num(i).eq.63) then                 
          v0=2*tt-s +2*hp+p
          u=2*xi-2*nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------Eta2-------------------------------------------------
        if(num(i).eq.64) then                 
          v0=2*tt+s+2*hp-p
          u=-2*nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------2Q1--------------------------------------------------
        if(num(i).eq.65) then                 
          v0=tt-4*s+hp+2*p+90.
          u=2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------Sigma1-------------------------------------------------
        if(num(i).eq.66) then                 
          v0=tt-4*s+3*hp+90.
          u=2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------Ro1-------------------------------------------------
        if(num(i).eq.67) then                 
          v0=tt-3*s+3*hp-p+90.
          u=2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------M11--------------------------------------------------
        if(num(i).eq.68) then                 
          v0=tt-s+hp+p-90.
          u=-nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------M12--------------------------------------------------
        if(num(i).eq.69) then                 
          v0=tt-s+hp-p-90.
          u=2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------Ki1--------------------------------------------------
        if(num(i).eq.70) then                 
          v0=tt-s+3*hp-p-90.
          u=-nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------Pi1--------------------------------------------------
        if(num(i).eq.71) then                 
          v0=tt-2*hp+p1+90.
          u=0.
          v0_u(i)=v0+u
          goto 130
        endif

c------------------Phi1-------------------------------------------------
        if(num(i).eq.72) then                 
          v0=tt+3*hp-90.
          u=0.
          v0_u(i)=v0+u
          goto 130
        endif

c------------------Teta1-------------------------------------------------
        if(num(i).eq.73) then                 
          v0=tt+s-hp+p-90.
          u=-nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------J1---------------------------------------------------
        if(num(i).eq.74) then                 
          v0=tt+s+hp-p-90.
          u=-nu
          v0_u(i)=v0+u
          goto 130
        endif

c------------------OO1--------------------------------------------------
        if(num(i).eq.75) then                 
          v0=tt+2*s+hp-90.
          u=-2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif
c
 130    v0_u(i)=amod(v0_u(i),360.00)    
c        write(*,*) mod(v0,360.),u
        

 140  continue

      return
      end


c########################################################################
c      
      subroutine nodalc(tj)
c
c------------------------------------------------------------------------
c  ROUTINE : nodalc
c
c  DESCRIPTION : Nodal correction calls
c
c
c  PROGRAMMERS : J.M. MOLINES and F. LYARD
c 
c  DATE : 19/01/94
c------------------------------------------------------------------------
      
      real*8 tj
      
      call astronomics(tj) 
      call nodal_a
      call nodal_G  

      return
      end


c########################################################################
c      
      subroutine init_data
c
c------------------------------------------------------------------------
c  ROUTINE : init_data
c
c  DESCRIPTION : Read the data files for each tidal constituents. The
C                path to the directory where the data files are located 
c                is given in wave_path variable
c                (e.g :wave_path='/users/toto/Grenoble_DATA').
c                nbext is the number of available constituents
c                Formerly 13 (FES95) Actually 27.
c                All files are named {wave_name}.model (ASCII files) or
c                {wave_name}.model.bin (BINARY files).
c                and model is a valid extension for the files
c
c  NOTE : August '2000 : No more admittance calculation for
c                 secondary constituents.
c
c  PROGRAMMERS : J.M. MOLINES - F. LYARD - F. LEFEVRE
c 
c  DATE : 15/06/1995 - 08/08/2000
c------------------------------------------------------------------------
      include 'common.h'  

      character*255 pathname, pathnameload
      character*20  name(29),extended(extended_max)
      real*8        frequency(29)
      integer*4     code(29),unit2 ,lex,lmod
      integer*4       lwa, strlen, nbext, nb, nn, istat, i

c-- Variables for admittance computation
      real*4 alk2,aln2,alm2,alnu2,almu2,all2,allda2,alt2 
      real*4 alq1,alo1,alk1,alp1
      integer*4 iq1,io1,ik1,im2,ik2,in2
      integer*4 inu2,imu2,il2,it2,ilda2,ip1
      real*4 frbar1,deno1,ck,sk,cn,sn
      real*4 deno,cnu2,snu2,cmu2,smu2,cl2,sl2,ct2,st2,clda2,slda2
      real*4 aap1,bbp1,ccp1
      real*4 aamu2,aanu2,aal2,aat2,aalda2,bbmu2,bbnu2,bbl2,bbt2,bblda2
      real*4 ccmu2,ccnu2,ccl2,cct2,cclda2
      common /admitancesaa/aamu2,aanu2,aal2,aat2,aalda2,aap1
      common /admitancesbb/bbmu2,bbnu2,bbl2,bbt2,bblda2,bbp1
      common /admitancescc/ccmu2,ccnu2,ccl2,cct2,cclda2,ccp1
c-- End of Variables for admittance computation

c
c these values are internal to the program and to the nodal correction
c routines. Do not change under any circumstances.
      data (code(i),name(i),frequency(i),i=1,8) /
     &     27,'Q1',13.39866087990d0,
     &     1,'O1' ,13.94303558000d0,
     &     3,'K1' ,15.04106864000d0,
     &     5,'2N2',27.89535481990d0,
     &     7, 'N2',28.43972952010d0,
     &     9, 'M2',28.98410422000d0,
     &     14,'K2',30.08213728000d0,
     &     13,'S2',30.00000000000d0/
c
      data (code(i),name(i),frequency(i),i=9,13) /
     &     2,'P1' ,14.95893136000d0,
     &     8,'Nu2',28.51258314000d0,
     &     6,'Mu2',27.96820844000d0,
     &     11,'L2',29.52847892000d0,
     &     12,'T2',29.95893332010d0/
c extra diurnal and semi diurnal are infered by admittance;
      data (code(i),name(i),frequency(i),i=14,28) /
     &     60,'Eps2',27.4238337d0,
     &     61,'Lambda2',29.4556253d0,
     &     64,'Eta2',30.6265120d0,
c
     &     65,'2Q1 ',12.8542862d0,
     &     66,'Sigma1',12.9271398d0,
     &     67,'Ro1',13.4715145d0,
     &     68,'M11 ',14.4966939d0,
     &     69,'M12 ',14.4966939d0,
     &     70,'Ki1 ',14.5695476d0,
     &     71,'Pi1 ',14.9178647d0,
     &     72,'Phi1',15.1232059d0,
     &     73,'Teta1',15.5125897d0,
     &     74,'J1  ',15.5854433d0,
     &     75,'OO1 ',16.1391017d0,
     &     76,'Psi1 ',15.0821353d0/

      lwa=strlen(wave_path)

      print *,'Initialisation ...please be patient ...'
      print *,'Used path variable :'
      print *,'path=',wave_path(1:lwa)
      if(iascbin.eq.1) then 
        print *,'Reading ASCII files'
      else
        print *,'Reading BIN files'
      endif

c
c
c Coefficient of the tidal potential (used in spline admittances)
c   spline admittances (see GRL 18(5):845-848,1991)
c
c
      alk2   = 0.1149327
      aln2   = 0.1758941
      alm2   = 0.9085024
      alnu2  = 0.03303
      almu2  = 0.02777
      all2   = 0.0251
      allda2 = 0.0066
      alt2   = 0.0247766
c
      alq1  = 0.073017
      alo1  = 0.3771366
      alk1  = 0.5300728
      alp1  = 0.1750754
c        
c internal index of the constituents
      iq1 = 1
      io1 = 2
      ik1 = 3
      im2 = 6
      ik2 = 7
      in2 = 5

      inu2  = 10
      imu2  = 11
      il2   = 12
      it2   = 13
      ilda2 = 15
      ip1   = 9

c compute the coefficient for P1 as a linear admittance with Q1,O1 K1
c (linear regression)
      frbar1=1/3.*(frequency(iq1) + 
     |     frequency(io1) + frequency(ik1))

      deno1=frbar1**2 -
     |     1/3.*(frequency(iq1)**2 + 
     |     frequency(io1)**2 + frequency(ik1)**2)
      aap1 = alp1/3./alq1*(1.-
     |     (frequency(ip1)-frbar1)*(frequency(iq1)-frbar1)/deno1)

      bbp1 = alp1/3./alo1*(1.-
     |     (frequency(ip1)-frbar1)*(frequency(io1)-frbar1)/deno1)

      ccp1 = alp1/3./alk1*(1.-
     |     (frequency(ip1)-frbar1)*(frequency(ik1)-frbar1)/deno1)

      ck=cos(2*pi*2.*(frequency(ik2)-frequency(im2))/15.) ! CPD
      sk=sin(2*pi*2.*(frequency(ik2)-frequency(im2))/15.) ! CPD

      cn=cos(2*pi*2.*(frequency(in2)-frequency(im2))/15.) ! CPD
      sn=sin(2*pi*2.*(frequency(in2)-frequency(im2))/15.) ! CPD
      deno=sk*(cn-1)-sn*(ck-1)
      
      cnu2=cos(2*pi*2.*(frequency(inu2)-frequency(im2))/15.)
      snu2=sin(2*pi*2.*(frequency(inu2)-frequency(im2))/15.)

      cmu2=cos(2*pi*2.*(frequency(imu2)-frequency(im2))/15.)
      smu2=sin(2*pi*2.*(frequency(imu2)-frequency(im2))/15.)

      cl2=cos(2*pi*2.*(frequency(il2)-frequency(im2))/15.)
      sl2=sin(2*pi*2.*(frequency(il2)-frequency(im2))/15.)

      ct2=cos(2*pi*2.*(frequency(it2)-frequency(im2))/15.)
      st2=sin(2*pi*2.*(frequency(it2)-frequency(im2))/15.)

      clda2=cos(2*pi*2.*(frequency(ilda2)-frequency(im2))/15.)
      slda2=sin(2*pi*2.*(frequency(ilda2)-frequency(im2))/15.)

      aamu2= (-sn*cmu2 +(cn-1)*smu2 +sn)/deno/alk2*almu2
      aanu2= (-sn*cnu2 +(cn-1)*snu2 +sn)/deno/alk2*alnu2
      aal2=  (-sn*cl2  +(cn-1)*sl2  +sn)/deno/alk2*all2
      aat2=  (-sn*ct2  +(cn-1)*st2  +sn)/deno/alk2*alt2
      aalda2=(-sn*clda2+(cn-1)*slda2+sn)/deno/alk2*allda2
c
      bbmu2=(sk*cmu2-(ck-1)*smu2-sk)/deno/aln2*almu2
      bbnu2=(sk*cnu2-(ck-1)*snu2-sk)/deno/aln2*alnu2
      bbl2=(sk*cl2-(ck-1)*sl2-sk)/deno/aln2*all2
      bbt2=(sk*ct2-(ck-1)*st2-sk)/deno/aln2*alt2
      bblda2=(sk*clda2-(ck-1)*slda2-sk)/deno/aln2*allda2
c
      ccmu2=(-(sk-sn)*cmu2+(ck-cn)*smu2+sk*cn-sn*ck)/deno/alm2*almu2
      ccnu2=(-(sk-sn)*cnu2+(ck-cn)*snu2+sk*cn-sn*ck)/deno/alm2*alnu2
      ccl2=(-(sk-sn)*cl2+(ck-cn)*sl2+sk*cn-sn*ck)/deno/alm2*all2
      cct2=(-(sk-sn)*ct2+(ck-cn)*st2+sk*cn-sn*ck)/deno/alm2*alt2
      cclda2=(-(sk-sn)*clda2+(ck-cn)*slda2+sk*cn-sn*ck)/deno/alm2*allda2
c
      nbext=8
c
      extended(1)='Q1'               
      extended(2)='O1'               
      extended(3)='K1'               
      extended(4)='2N2' 
      extended(5)='N2'               
      extended(6)='M2'               
      extended(7)='K2'               
      extended(8)='S2' 
      extended(9)='P1'               
      extended(10)='Nu2'               
      extended(11)='Mu2'               
      extended(12)='L2'               
      extended(13)='T2'
      extended(14)='Eps2'
      extended(15)='Lambda2'
      extended(16)='Eta2'
      extended(17)='2Q1'
      extended(18)='Sigma1'
      extended(19)='Ro1'
      extended(20)='M11'
      extended(21)='M12'
      extended(22)='Ki1'
      extended(23)='Pi1'
      extended(24)='Phi1'
      extended(25)='Teta1'
      extended(26)='J1'
      extended(27)='OO1'
      extended(28)='Psi1'
      
      unit2=20
      
      nb=0      
 742  continue
      do nn=1,nbext 

        nb=nb+1 
        lex=strlen(extended(nn))
        lmod=strlen(model)
      if(iascbin.eq.1) then
        pathname=wave_path(1:lwa)//
     &       extended(nn)(1:lex)//'_'//model(1:lmod)//'.asc'
        pathnameload=wave_path(1:lwa)//
     &       extended(nn)(1:lex)//'_dr'//model(1:lmod)//'.asc'

      else
        pathname=wave_path(1:lwa)//
     &       extended(nn)(1:lex)//'_'//model(1:lmod)//'.bin'
        pathnameload=wave_path(1:lwa)//
     &       extended(nn)(1:lex)//'_dr'//model(1:lmod)//'.bin'
      endif
c
c
 120    wave(nb)=extended(nn)
        print *,'Read file : ',pathname(1:strlen(pathname))
        print *,'Read file : ',pathnameload(1:strlen(pathnameload))

        call read_tide(nb,unit2,pathname,pathnameload,istat)

        if (istat.eq.0) then 
          nbwave=nb
          do i=1,9
            if(wave(nb).eq.name(i)) then
              freq(nb)=frequency(i)
              num(nb)=code(i)   
              freq(nb)=freq(nb)*pi/180.d+00
              goto 100
            endif
          enddo 
          stop 'Unknown wave...'
 100      continue
        else
          nb=nb-1
        endif
 101    continue
      enddo  

 10   format(i3,3x,a6,3x,d14.11) 
c
c nine (8) wave are read from the files 19 are inferred by admittance
c
      nbwave=27
      do i=9,nbwave   
        nb=nb+1   
        wave(nb)=extended(i)
        freq(nb)=frequency(i)
        num(nb)=code(i)   
        freq(nb)=freq(nb)*pi/180.d+00
      enddo 
      
      call init_corrections(0.0d0)
      print *, 'Now ready !!!'
      return
      end


c########################################################################
c      
      subroutine read_tide(nb,iunit,pathname,pathnameload,istat)
c
c------------------------------------------------------------------------
c  ROUTINE :read_tide
c
c  DESCRIPTION : Read data files either  ASCII formatted or BIN files
c                (Check the CPPFLAG in Makefile)
c                to convert ASCII to BIN use ascii2bin.f program
c
c  PROGRAMMERS :  F. LYARD and F. LEFEVRE
c 
c  DATE : 11/11/2000
c------------------------------------------------------------------------
      include 'common.h'   

      integer*4 iunit
      integer*4 nb, istat, i, j
      character*255 pathname, pathnameload
      integer*4 nimax,njmax
      real*4 conv
      real*4 xmax, ymax
      real*4 amask
      real*4 Gmask
      integer*4 k
      parameter (nimax=1441,njmax=721)
      real*4 wra(nimax,njmax), wrg(nimax,njmax), dn
      character*80 comment(4)
      integer*4 nk, nt, nd, icode
      real*4 mask       
      real*4 a,G
      real*4 level
      
c
c default files are ASCII files, (BIN_FILE not defined)
c
      if(iascbin.eq.1) then    

c### Read ASCII file format

        conv=pi/180. 
        spec=32767.
        istat=0

c--- Pure oceanic tide
        open(iunit,file=pathname,status='old',err=100)
        read(iunit,*) xmin,xmax
        read(iunit,*) ymin,ymax
        read(iunit,*) dx,dy
        read(iunit,*) ni,nj
        read(iunit,*) amask,Gmask

        do j=1,nj
          do k=1,ni,30
            read(iunit,*) (wra(i,j),i=k,MIN(ni,k+29))
            read(iunit,*) (wrg(i,j),i=k,MIN(ni,k+29))
          enddo
        enddo

        do i=1,ni
          do j=1,nj
            if (wra(i,j).ne.amask.and.wrg(i,j).ne.Gmask) then
              wrp(i,j,nb)=wra(i,j)*cos(-wrg(i,j)*conv)
              wip(i,j,nb)=-wra(i,j)*sin(-wrg(i,j)*conv)
            else
              wrp(i,j,nb)=spec
              wip(i,j,nb)=spec
            endif
          enddo

        enddo          

        close(iunit) 

c--- Radial loading tide
        open(iunit,file=pathnameload,status='old',err=100)
        read(iunit,*) xmin,xmax
        read(iunit,*) ymin,ymax
        read(iunit,*) dx,dy
        read(iunit,*) ni,nj
        read(iunit,*) amask,Gmask

        do j=1,nj
          do k=1,ni,30
            read(iunit,*) (wra(i,j),i=k,MIN(ni,k+29))
            read(iunit,*) (wrg(i,j),i=k,MIN(ni,k+29))
          enddo
        enddo
 
c--- Shift 0 to 360 degrees
        dn=-xmin/dx
        do i=1,ni
          if(i.le.dn) then
            k=i+dn+1
          else
            k=i-dn
          endif
          do j=1,nj
            if (wra(i,j).ne.amask.and.wrg(i,j).ne.Gmask) then
              wrpload(k,j,nb)=wra(i,j)*cos(wrg(i,j)*conv)
              wipload(k,j,nb)=wra(i,j)*sin(wrg(i,j)*conv)
            else
              wrpload(k,j,nb)=spec
              wipload(k,j,nb)=spec
            endif
          enddo

        enddo          

c--- Shift 0 to 360 degrees
        xmin=0.0

        close(iunit) 

      else   

c### Read BIN file format : compatibility with LEGOS Bimg format

        spec=32767.
        istat=0

c--- Pure oceanic tide
        open(iunit,file=pathname,status='old',
     &       form='unformatted',err=100)

        do i=1,4
          read(iunit) comment(i)
        enddo
        read(iunit) ni,nj,nk,nt,nd
        read(iunit) xmin,ymin,dx,dy,mask,icode
        read(iunit) level
        read(iunit) level

        read(iunit) ((wrp(i,j,nb),wip(i,j,nb),i=1,ni),j=1,nj)

        close(iunit)

c--- Radial loading tide
        open(iunit,file=pathnameload,status='old',
     &       form='unformatted',err=100)

        do i=1,4
          read(iunit) comment(i)
        enddo
        read(iunit) ni,nj,nk,nt,nd
        read(iunit) xmin,ymin,dx,dy,mask,icode
        read(iunit) level
        read(iunit) level

        read(iunit) ((wrpload(i,j,nb),wipload(i,j,nb),i=1,ni),j=1,nj)

        close(iunit)

c Convert in cm => x100.
        do j=1,nj
          do i=1,ni
            a=wrp(i,j,nb)
            G=wip(i,j,nb)
            if (a.ne.mask.and.G.ne.mask) then
              wrp(i,j,nb)=wrp(i,j,nb)*100.
              wip(i,j,nb)=-wip(i,j,nb)*100.
              wrpload(i,j,nb)=wrpload(i,j,nb)*100.
              wipload(i,j,nb)=-wipload(i,j,nb)*100.
            else
              wrp(i,j,nb)=spec
              wip(i,j,nb)=spec
              wrpload(i,j,nb)=spec
              wipload(i,j,nb)=spec
            endif
          enddo
        enddo      


      endif
      
      print *,'  Wave ',wave(nb),' ok...' 
      return
      
 100  continue
      print *,' ERROR opening data file', pathname
      istat=-1
      stop 'FATAL ERROR'

      end


c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c    FUNCTIONS
c    
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



c########################################################################
c      
      real*8  function   deg(z)   
c
c------------------------------------------------------------------------
c  ROUTINE : deg
c
c  DESCRIPTION : conversion from radians to deg
c
c  PROGRAMMERS : F. LYARD
c 
c  DATE : 19/01/94
c------------------------------------------------------------------------
      include 'common.h'          
      real*8 z    
                          
      deg=(z*1.8D+02)/pi
      return
      end


c########################################################################
c      
      real*8 function rad(w)      
c
c------------------------------------------------------------------------
c  ROUTINE : rad
c
c  DESCRIPTION : conversion de degres en radians
c
c  PROGRAMMERS : J.M. MOLINES et F. LYARD
c 
c  DATE : 19/01/94
c------------------------------------------------------------------------
       include 'common.h'          
c------------------------------------------------------------------------
c
c     cette fonction permet de convertir la mesure d'un angle exprimee
c     en degres --> radians.
c                         
c------------------------------------------------------------------------
      real*8 w    
      rad=(w*pi)/1.8D+02

      return
      end


c########################################################################
c      
      real*8 function julian(annee,jour,seconde)
c------------------------------------------------------------------------
c  ROUTINE :julian
c
c  DESCRIPTION : This function return the elapsed time 
c                since Jan 1 1900 0:0TU
c                in julian centuries .
c
c  PROGRAMMERS : J.M. MOLINES et F. LYARD
c 
c  DATE : 19/01/94
c------------------------------------------------------------------------
c
      real*8  xj,bissex
      integer*4 annee,jour
      real*4    seconde

      xj=dble((annee-1900.)*365)
      xj=xj+dble(jour-1)
      xj=xj+dble(seconde)/8.64d+04
      bissex=dble(int((annee-1901)/4))
      julian=(xj+bissex)/3.6525d+04  
      
      return
      end



c########################################################################
c      
      function strlen(string)
c
c------------------------------------------------------------------------
c  ROUTINE : strlen
c
c  DESCRIPTION : This function is usually part of libU77 on unix system. 
c                provided for portability
c
c  PROGRAMMERS : J.M. MOLINES
c 
c  DATE : 19/01/94
c------------------------------------------------------------------------
      character string*(*)
      integer*4 ll, ii, strlen, k

      ll= len(string)
      ii= index(string,' ')
      if (ll.eq.0.or.ii.eq.1) then
        strlen= 1
        return
      else if (ii.eq.0) then
        strlen= ll
        return
      end if
      do 10 k= ll, ii-1, -1
        if (string(k:k).ne.' ') go to 20
10        continue
20      strlen= k
      return
      end


c########################################################################
c      
      integer*4 function tojul (ian,imois,ijour)
c
c------------------------------------------------------------------------
c  ROUTINE : tojul
c
c  DESCRIPTION : renvoi une date au format aajjj quand elle est donnee
c                au format aaaa mm jj 
c
c  PROGRAMMERS : J.M. MOLINES and F. LEFEVRE
c 
c  DATE  : 19/01/94 - 08/08/2000
c------------------------------------------------------------------------
      dimension   jmois(12)
      integer*4   ian, imois, ijour, jmois, ia, m, ida
      real*4 s                    
      data jmois /31,28,31,30,31,30,31,31,30,31,30,31/
      
      if(mod(ian,4).eq.0) then
        jmois(2)=29
      else
        jmois(2)=28
      endif
      if(ijour.gt.jmois(imois)) then
        print *,' ERROR in tojul : There are only ',jmois(imois),
     &       'days in month number ',imois
        stop 'FATAL ERROR'
      endif
      ia=ian-1900
      s=0
      do m=1,imois-1
        s=s+jmois(m)
      enddo
      ida=s+ijour
      tojul=ia*1000+ida
      return
      end


