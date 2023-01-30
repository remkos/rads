#ifndef GRID_FILE
#ifndef BIMG_FILE
#define ASCII_FILE
#endif
#endif
c###
	subroutine otide(model_name,tide,lat,long,jnasa,heure,istat)
c------------------------------------------------------------------------------
c  PROGRAMME : otide.f
c
c  DESCRIPTION : tide prediction routine (end user routine):
c	   tide= predicted tide in cm (output)
c          lat = latitude (between -90 et +90 deg ) (input)
c          long = longitude (either between -180 180 or 0 360).
c          jnasa= Day for prediction (as in TP GDR )
c          heure= Hour during the day jnasa
c          istat= Number of valid points in the interpolation
c
c  PROGRAMMEURS : J.M. MOLINES et Remko SCHARROO (DUT/SSR&T)
c 
c  DATE derniere retouche : 16/06/95
c  Majeur revision : 11 Jul 1995 par Remko Scharroo.
c------------------------------------------------------------------------------
        include 'common.h' 
c
	real*8 t1,heure
	character*(*) model_name
	real*4 lat,long,tide,lon
	logical first
	data first /.true./
	save first

c read the data files at first call. This is quite long
	if (first) then
	model=model_name
	 call init_data
	 first=.false.
	endif
c
c Conversion  TOPEX date to adequate date for prediction
c * TOPEX date is MJD since 1.0 Jan 1958
c * T1 is date in Julian Centuries since 1.0 Jan 1900
c * 21184 = days from 1900 to 1958
c * 36525 = number of days in a Julian Century
c
        t1=(jnasa+21184+heure/24d0)/36525d0
c
c Convert longitude in range -540 to +540 to range -180 to +180
c
	if (long.ge.180.) then
	   lon=long-360.
	else if (long.lt.-180.) then
	   lon=long+360.
	else
	   lon=long
	endif

        call grenoble(t1,lon,lat,tide,istat)

	return
	end
c##
                subroutine astronomics(tj)                                          
c------------------------------------------------------------------------------
c  PROGRAMME : astronomics.f
c
c  DESCRIPTION : This program initialize some astronomic data usefull for
c                nodal corrections.
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
       include 'common.h'          
                      
      real*8 tj,tgn2,at1,at2,u,tgi2
      real*8 ct0,ct1
      real*8 cn0,cn1
      real*8 cs0,cs1
      real*8 ch0,ch1
      real*8 cps0,cps1
      real*8 cp0,cp1  

c------------- tt mean solar angle relative to Greenwich ---------------
                    
	ct0=180.D+00
	ct1=360.D+00*3.6525D+04
      tt=dmod(ct0+ct1*tj,360.D+00)

c-------------- hp longitude of ascending lunar node ----------------
                                                           
	cn0=  259.1560563D+00
	cn1= -1934.1423972D+00   
      n=dmod(cn0+cn1*tj,360.D+00)

c------------------- hp mean solar longitude ----------------------

	ch0= 280.1895015D+00
	ch1= 36000.76892D+00 
	hp=dmod(ch0+ch1*tj,360.D+00)
                           
c-------------------- s mean lunar longitude ----------------------

	cs0= 277.0256206D+00
	cs1= 481267.892D+00
	s=dmod(cs0+cs1*tj,360.D+00)

c------------------- p1 longitude of solar perigee --------------------
                     
	cps0=281.2208568D+00
	cps1=1.719175D+00   
	p1=dmod(cps0+cps1*tj,360.D+00)

c------------------- p longitude of lunar perigee --------------------

	cp0=334.3837214D+00
	cp1=4069.0322056D+00
	p=dmod(cp0+cp1*tj,360.D+00)

      n=n*rad
      hp=hp*rad
      s=s*rad
      p1=p1*rad
      p=p*rad

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
   
c------------------------for constituents l2,k1,k2------------------------


      tgi2=dtan(iang/2.D+00)     
      pp=p-xi
      x1ra=dsqrt(1.-12.D+00*tgi2**2*cos(2.D+00*pp)+36.D+00*tgi2**4)
      r=datan(dsin(2.*pp)/(1./(6.*tgi2**2)-dcos(2.*pp)))
      nuprim=datan(dsin(2.D+00*iang)*dsin(nu)/
     &            (dsin(2.D+00*iang)*dcos(nu)+3.347D-01)) 
      nusec=0.5*datan(((dsin(iang)**2.D+00)*dsin(2.D+00*nu))
     &   /(dsin(iang)**2.D+00*dcos(2.*nu)+7.27D-02))                 

c------------------------------------------------------------------------
      return
      end
c###
      subroutine grenoble(t,rlon,rlat,vtide,istat)
c------------------------------------------------------------------------------
c  PROGRAMME : grenoble.f
c
c  DESCRIPTION : driver for tidal computation. The result is returned in tide.
c               tide is an array of 
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 13/03/95
c  Bug removed (caused by 'newpi') : 19/02/97 : Remko SCHARROO
c------------------------------------------------------------------------------
       include 'common.h' 
       dimension tide(extended_max)        

c----------------------------------------------------------------------- 
c istat = Number of grid points for data interpolation
c istat = 0  point (rlon,rlat) is out of the gridded area
c istat = 1,2,3  point (rlon,rlat) is on a frontier of the area
c istat = 4  point (rlon,rlat) is fully in the domain
c-----------------------------------------------------------------------

	real*8 delta,t
                  
      delta=(t-t_nodal)*3.6525d+04*24.D+00

	if(abs(delta).gt.delta_max) then
	  call init_corrections(t)
	  delta=0.D+00
	endif

c-----------------------------------------------------------------------
c      sreal(nb) = real part of the nb th tide wave
c      simag(nb) = imaginary  part of the nb th tide wave
c-----------------------------------------------------------------------

	call interpolation(rlon,rlat,istat)
                     
	if (istat.ne.0) then
	  h=0.0
	 
         do nb=1,nbwave
          phi=real(freq(nb)*delta)+v0_u(nb)  
          h=h+f(nb)*(sreal(nb)*cos(phi)+simag(nb)*sin(phi))
	  tide(nb)=h
c        print 1956,nb,wave(nb),sreal(nb),simag(nb),sqrt(sreal(nb)**2+
c    |      simag(nb)**2),v0_u(nb),phi,freq(nb),delta
c1956	format(i3,1x,a10,7f8.3)
c23456
         enddo 
	else
	 do nb=1,nbwave
	  tide(nb)=spec
	 enddo
	endif
	vtide=tide(nbwave)
    
      return
      end
c###
      subroutine init_corrections(t0)
c------------------------------------------------------------------------------
c  PROGRAMME : init_corrections.f
c
c  DESCRIPTION : compute nodal corrections
c
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
       include 'common.h'          
                           
      real*8 t0
                  
      call nodalc(t0)     
                   
      do nb=1,nbwave   
c-->    convert v0_u from degrees to radians
	  v0_u(nb)=v0_u(nb)*rad
      enddo   

	t_nodal=t0

      return
      end
c###
      subroutine interpolation(xx,yy,istat)
c------------------------------------------------------------------------------
c  PROGRAMME : interpolation.f
c
c  DESCRIPTION : Perform bi-linear interpolation at point xx,yy from the gridded
c                files.
c          istat returns the number of points used for the interpolation
c   June 1995: new constituents added by admittance after interpolation
c              for the main constituents
c
c
c  PROGRAMMEURS : J.M. MOLINES 
c 
c  DATE derniere retouche : 15/06/1995
c------------------------------------------------------------------------------
       include 'common.h'   

	common /admitancesaa/aamu2, aanu2, aal2, aat2,aalda2,aap1
	common /admitancesbb/bbmu2, bbnu2, bbl2, bbt2,bblda2,bbp1
	common /admitancescc/ccmu2, ccnu2, ccl2, cct2,cclda2,ccp1
                       
	istat=0     

	if (nbwave.eq.0) return

      do nb=1,nbwave 
        sreal(nb)=0.0
        simag(nb)=0.0  
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
	    if(wrp(i,j,nb).gt.(spec-0.1).
     .        or.wip(i,j,nb).gt.(spec-0.1)) goto 200
	    end do
c all minor constituents will be deduced from the major ones.
c
	    istat=istat+1
	    xij=xmin+(i-1)*dx  
	    yij=ymin+(j-1)*dy
c this line modified for taking into account interpolation near the origin
c           pds=(dx-amod(abs(x-xij),dx))*(dy-abs(y-yij))
            pds=(dx-abs(x-xij))*(dy-abs(y-yij))
	    ptot=ptot+pds   
          do nb=1,8
            sreal(nb)=sreal(nb)+wrp(i,j,nb)*pds    
            simag(nb)=simag(nb)+wip(i,j,nb)*pds
	  enddo
200	    continue
        enddo 
100	  continue 
      enddo

	do nb=1,nbwave
	  if(ptot.ne.0.0) then
	    sreal(nb)=sreal(nb)/ptot  
	    simag(nb)=simag(nb)/ptot  
	  else    
	    istat=0
	    sreal(nb)=spec  
	    simag(nb)=spec  
	  endif
        enddo 
c
	if (istat.ne.0) then
c
c infer additional constituents by admittance
c
c DIURNALS (from Richard Ray perth2 program. Thank You, Richard!)
c---------
c  from Q1 and O1 (1-2)
c    2Q1
          sreal(17) = 0.263 *sreal(1) - 0.0252*sreal(2)
          simag(17) = 0.263 *simag(1) - 0.0252*simag(2)
c    sigma1
          sreal(18) = 0.297 *sreal(1) - 0.0264*sreal(2)
          simag(18) = 0.297 *simag(1) - 0.0264*simag(2)
c    rho1
          sreal(19) = 0.164 *sreal(1) + 0.0048*sreal(2)
          simag(19) = 0.164 *simag(1) + 0.0048*simag(2)
c  from O1 and K1  (2-3)
c    M11
          sreal(20) = 0.0140*sreal(2) + 0.0101*sreal(3)
          simag(20) = 0.0140*simag(2) + 0.0101*simag(3)
c    M12
          sreal(21) = 0.0389*sreal(2) + 0.0282*sreal(3)
          simag(21) = 0.0389*simag(2) + 0.0282*simag(3)
c    chi1
          sreal(22) = 0.0064*sreal(2) + 0.0060*sreal(3)
          simag(22) = 0.0064*simag(2) + 0.0060*simag(3)
c    pi1
          sreal(23) = 0.0030*sreal(2) + 0.0171*sreal(3)
          simag(23) = 0.0030*simag(2) + 0.0171*simag(3)
c    phi1
          sreal(24) =-0.0015*sreal(2) + 0.0152*sreal(3)
          simag(24) =-0.0015*simag(2) + 0.0152*simag(3)
c    theta1
          sreal(25) =-0.0065*sreal(2) + 0.0155*sreal(3)
          simag(25) =-0.0065*simag(2) + 0.0155*simag(3)
c   J1
          sreal(26) =-0.0389*sreal(2) + 0.0836*sreal(3)
          simag(26) =-0.0389*simag(2) + 0.0836*simag(3)
c   OO1
          sreal(27) =-0.0431*sreal(2) + 0.0613*sreal(3)
          simag(27) =-0.0431*simag(2) + 0.0613*simag(3)
c
c    P1  from Grenoble admittance code
	  sreal(9) =  aap1*sreal(1)+bbp1*sreal(2)+ccp1*sreal(3)
	  simag(9) =  aap1*simag(1)+bbp1*simag(2)+ccp1*simag(3)
 
c SEMI-DIURNAL (from Grenoble to take advantage of 2N2 )
c ------------
c  from 2N2 -N2 (4-5)
c    eps2 
          sreal(14) = 0.53285 *sreal(4) - 0.03304*sreal(5)
          simag(14) = 0.53285 *simag(4) - 0.03304*simag(5)
c  from M2 - K2 (6-7)
c   eta2
          sreal(16) = -0.0034925 *sreal(6) + 0.0831707*sreal(7)
          simag(16) = -0.0034925 *simag(6) + 0.0831707*simag(7)
c
c  from N2 -M2- K2 by spline admittances (see GRL 18(5):845-848,1991)
c
c   mu2
          sreal(11) = aamu2*sreal(7)+bbmu2*sreal(5)+ccmu2*sreal(6)
          simag(11) = aamu2*simag(7)+bbmu2*simag(5)+ccmu2*simag(6)
c   nu2
          sreal(10) = aanu2*sreal(7)+bbnu2*sreal(5)+ccnu2*sreal(6)
          simag(10) = aanu2*simag(7)+bbnu2*simag(5)+ccnu2*simag(6)
c   lda2
          sreal(15) = aalda2*sreal(7)+bblda2*sreal(5)+cclda2*sreal(6)
          simag(15) = aalda2*simag(7)+bblda2*simag(5)+cclda2*simag(6)
c   L2
          sreal(12) = aal2*sreal(7)+bbl2*sreal(5)+ccl2*sreal(6)
          simag(12) = aal2*simag(7)+bbl2*simag(5)+ccl2*simag(6)
c   T2
          sreal(13) = aat2*sreal(7)+bbt2*sreal(5)+cct2*sreal(6)
          simag(13) = aat2*simag(7)+bbt2*simag(5)+cct2*simag(6)

	endif
    
      return
      end
*c###
*	function lnblnk(string)
*c------------------------------------------------------------------------------
*c  PROGRAMME : lnblnk.f
*c
*c  DESCRIPTION : This function is usually part of libU77 on unix system. 
*c                 provided for portability
*c
*c  PROGRAMMEURS : J.M. MOLINES
*c 
*c  DATE derniere retouche : 19/01/94
*c------------------------------------------------------------------------------
*	character string*(*)
*	ll= len(string)
*	ii= index(string,' ')
*	if (ll.eq.0.or.ii.eq.1) then
*	  lnblnk= 1
*	  return
*	else if (ii.eq.0) then
*	  lnblnk= ll
*	  return
*	end if
*	do 10 k= ll, ii-1, -1
*	  if (string(k:k).ne.' ') go to 20
*10	  continue
*20	lnblnk= k
*	return
*	end
*c###
         subroutine nodal_a      
c------------------------------------------------------------------------------
c  PROGRAMME : nodal_a.f
c
c  DESCRIPTION :Compute nodal corrections from SCHUREMAN (1958)
c     
c  CAUTION     : indexes used in this routine are internal to the code
c                and corresponds to the !original! ondes.dat file.
c      june 1995: new secondary constituents added.
c
c  PROGRAMMEURS : J.M. MOLINES
c 
c  DATE derniere retouche : 15/06/1995
c------------------------------------------------------------------------------
       include 'common.h'     
     
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
11        f(i)=dsin(iang)*dcos(iang/2)**2/0.38  
          goto 120
12        f(i)=1 
          goto 120 
13        f(i)=sqrt(0.8965*dsin(2*iang)**2+0.6001*dsin(2*iang)
     &              *dcos(nu)+0.1006)          
          goto 120
14        f(i)=real(dcos(iang/2.D+00)**4/9.154D-01)
          goto 120
15        f(i)=dcos(iang/2)**4/0.9154*x1ra        
          goto 120 
16        f(i)=sqrt(19.0444*dsin(iang)**4
     &               +2.7702*dsin(iang)**2*dcos(2*nu)+0.0981)
	  goto 120
17	  f(i)=real(dsin(iang)*dsin(iang)/0.1565d0) 
          goto 120      
18	  f(i)=real(dsin(2*iang)/0.7214d0) 
          goto 120
19	  f(i)=real(dsin(iang)*dsin(iang/2.d0)*dsin(iang/2.d0)/0.0164d0) 
          goto 120
      
120     continue   
100   continue
      return
      end
c###
                  subroutine  nodal_G  
c------------------------------------------------------------------------------
c  PROGRAMME : nodal_g.f
c
c  DESCRIPTION : compute V0+u from Schureman (1958)
c       june 1995: new secondary constituents added
c
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 16/06/95
c------------------------------------------------------------------------------
       include 'common.h'          

      n=n/rad
      hp=hp/rad
      s=s/rad
      p1=p1/rad
      p=p/rad
      xi=xi/rad
      nu=nu/rad
      nuprim=nuprim/rad
      nusec=nusec/rad
      r=r/rad

      do 140 i=1,nbwave
c-------------------o1-------------------------------------------------- 
        if(num(i).eq.1) then
          v0=tt-2*s+hp+90.
          u=2*xi-nu
          v0_u(i)=tt-2*s+hp+90.+2*xi-nu
          goto 130
        endif
c-------------------p1--------------------------------------------------
        if(num(i).eq.2) then  
          v0=tt-hp+90.
          u=0
          v0_u(i)=tt-hp+90.
          goto 130
        endif      
c-------------------k1--------------------------------------------------
        if(num(i).eq.3) then           
          v0=tt+hp-90.                                             
          u=-nuprim
          v0_u(i)=tt+hp-90.-nuprim
          goto 130 
        endif
                                   
c-------------------2n2------------------------------------------------
        if(num(i).eq.5) then
          v0=2*tt-4*s+2*hp+2*p
          u=2*xi-2*nu
          v0_u(i)=2*tt-4*s+2*hp+2*p+2*xi-2*nu
          goto 130
        endif                                                           
c-------------------mu2-------------------------------------------------
        if(num(i).eq.6) then
          v0=2*tt-4*s+4*hp
          u=2*xi-2*nu
          v0_u(i)=2*tt-4*s+4*hp+2*xi-2*nu
          goto 130 
        endif
c-------------------n2--------------------------------------------------
       if(num(i).eq.7) then
          v0=2*tt-3*s+2*hp+p
          u=2*xi-2*nu
          v0_u(i)=2*tt-3*s+2*hp+p+2*xi-2*nu
          goto 130
        endif
c-------------------nu2-------------------------------------------------
        if(num(i).eq.8) then
          v0=2*tt-3*s+4*hp-p
          u=2*xi-2*nu
          v0_u(i)=2*tt-3*s+4*hp-p+2*xi-2*nu
          goto 130
        endif     
c-------------------m2--------------------------------------------------
        if(num(i).eq.9) then
          v0=2.*tt-2*s+2*hp
          u=2*xi-2*nu
          v0_u(i)=2.*tt+real(-2.D+00*s+2.D+00*hp+2.D+00*xi-2.D+00*nu)
          goto 130
        endif
c------------------l2---------------------------------------------------
        if(num(i).eq.11) then
          v0=2*tt-s+2*hp-p+180
          u=2*xi-2*nu-r
          v0_u(i)=2*tt-s+2*hp-p+180+2*xi-2*nu-r
          goto 130
        endif
c------------------t2---------------------------------------------------
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
c------------------k2---------------------------------------------------
        if(num(i).eq.14) then
          v0=2*tt+2*hp
          u=-2*nusec
          v0_u(i)=2*tt+2*hp-2*nusec
          goto 130
        endif
c------------------q1---------------------------------------------------
        if(num(i).eq.27) then                 
          v0=tt-3*s+hp+p+90.
          u=2*xi-nu
          v0_u(i)=tt-3*s+hp+p+90.+2*xi-nu
          goto 130
        endif
c------------------eps2---------------------------------------------------
        if(num(i).eq.60) then                 
          v0=2*tt-5*s+4*hp+p
          u=0.
          v0_u(i)=v0+u
          goto 130
        endif
c------------------lda2---------------------------------------------------
        if(num(i).eq.61) then                 
          v0=2*tt-s+p+180.
          u=2*xi-2*nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------l21---------------------------------------------------
        if(num(i).eq.62) then                 
          v0=2*tt-s +2*hp-p+180
          u=-2*nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------l22---------------------------------------------------
        if(num(i).eq.63) then                 
          v0=2*tt-s +2*hp+p
          u=2*xi-2*nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------eta2---------------------------------------------------
        if(num(i).eq.64) then                 
          v0=2*tt+s+2*hp-p
          u=-2*nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------2q1---------------------------------------------------
        if(num(i).eq.65) then                 
          v0=tt-4*s+hp+2*p+90.
          u=2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------sig1---------------------------------------------------
        if(num(i).eq.66) then                 
          v0=tt-4*s+3*hp+90.
          u=2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------rho1---------------------------------------------------
        if(num(i).eq.67) then                 
          v0=tt-3*s+3*hp-p+90.
          u=2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------m11---------------------------------------------------
        if(num(i).eq.68) then                 
          v0=tt-s+hp+p-90.
          u=-nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------m12---------------------------------------------------
        if(num(i).eq.69) then                 
          v0=tt-s+hp-p-90.
          u=2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------xi1---------------------------------------------------
        if(num(i).eq.70) then                 
          v0=tt-s+3*hp-p-90.
          u=-nu
          v0_u(i)=v0+u
          goto 130
        endif
c------------------pi1---------------------------------------------------
        if(num(i).eq.71) then                 
          v0=tt-2*hp+p1+90.
          u=0.
          v0_u(i)=v0+u
          goto 130
        endif
c------------------phi1---------------------------------------------------
        if(num(i).eq.72) then                 
          v0=tt+3*hp-90.
          u=0.
          v0_u(i)=v0+u
          goto 130
        endif
c------------------tta1---------------------------------------------------
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
c------------------OO1---------------------------------------------------
        if(num(i).eq.75) then                 
          v0=tt+2*s+hp-90.
          u=-2*xi-nu
          v0_u(i)=v0+u
          goto 130
        endif
c
130     v0_u(i)=amod(v0_u(i),360.00)    
c	  write(*,*) mod(v0,360.),u
	  

140   continue

      return
      end
c###
               subroutine nodalc(tj)
c------------------------------------------------------------------------------
c  PROGRAMME : nodalc.f
c
c  DESCRIPTION : Nodal correction calls
c
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
                    
      real*8 tj
                                                                
      call astronomics(tj) 
      call nodal_a
      call nodal_G  

	return
      end
c###
      subroutine init_data
c------------------------------------------------------------------------------
c  PROGRAMME : init_data
c
c  DESCRIPTION : Read the data files for each tidal constituents. The path
c               to the directory where the data files are located is given
c               by an environment variable WAVE_PATH, that must be setenv
c               before running the program.
c              (e.g :setenv WAVE_PATH /users/toto/Grenoble_DATA ).
c               nbext is the number of available constituents. Actually 13.
c              All files are named {wave_name}.model (ASCII files) or
c                                  {wave_name}.model.bimg (BINARY files).
c                       and model is a valid extension for the files
c  June 1995 : Modified to prepare the admittance calculation for
c             secondary constituents.
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 15/06/1995
c------------------------------------------------------------------------------
       include 'common.h'  

      character*255    pathname    

      character*6      name(29),extended(extended_max)
      character*3 rep
      character*255 WAVE_PATH  
      real*8 frequency(29)
      integer*4        code(29),unit2 ,lex,lmod 

	common /admitancesaa/aamu2, aanu2, aal2, aat2,aalda2,aap1
	common /admitancesbb/bbmu2, bbnu2, bbl2, bbt2,bblda2,bbp1
	common /admitancescc/ccmu2, ccnu2, ccl2, cct2,cclda2,ccp1
c
c these values are internal to the program and to the nodal correction
c routines. Do not change under any circumstances.
	data (code(i),name(i),frequency(i),i=1,8) /
     | 27,'q1',13.39866087990d0,
     | 1,'o1' ,13.94303558000d0,
     | 3,'k1' ,15.04106864000d0,
     | 5,'2n2',27.89535481990d0,
     | 7, 'n2',28.43972952010d0,
     | 9, 'm2',28.98410422000d0,
     | 14,'k2',30.08213728000d0,
     | 13,'s2',30.00000000000d0/
c
	data (code(i),name(i),frequency(i),i=9,13) /
     | 2,'p1' ,14.95893136000d0,
     | 8,'nu2',28.51258314000d0,
     | 6,'mu2',27.96820844000d0,
     | 11,'l2',29.52847892000d0,
     | 12,'t2',29.95893332010d0/
c extra diurnal and semi diurnal are infered by admittance;
       data (code(i),name(i),frequency(i),i=14,27) /
     | 60,'eps2',27.4238337d0,
     | 61,'lda2',29.4556253d0,
     | 64,'eta2',30.6265120d0,
c
     | 65,'2q1 ',12.8542862d0,
     | 66,'sig1',12.9271398d0,
     | 67,'rho1',13.4715145d0,
     | 68,'m11 ',14.4966939d0,
     | 69,'m12 ',14.4966939d0,
     | 70,'xi1 ',14.5695476d0,
     | 71,'pi1 ',14.9178647d0,
     | 72,'phi1',15.1232059d0,
     | 73,'tta1',15.5125897d0,
     | 74,'j1  ',15.5854433d0,
     | 75,'oo1 ',16.1391017d0/
c
c Coefficient of the tidal potential (used in spline admittances)
c   spline admittances (see GRL 18(5):845-848,1991)
c
c
	alk2   = 0.1149327
	aln2   = 0.1758941
	alm2   = 0.9085024
	alnu2 = 0.03303
	almu2 = 0.02777
	all2  = 0.0251
	allda2= 0.0066
	alt2  = 0.0247766
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
	im2=6
	ik2=7
	in2=5

	inu2=10
	imu2=11
	il2=12
	it2=13
	ilda2=15
	ip1 = 9
c compute the coefficient for P1 as a linear admittance with Q1,O1 K1
c (linear regression)
	frbar1=1/3.*(frequency(iq1) + 
     |             frequency(io1) + frequency(ik1))

	deno1=frbar1**2 -
     |        1/3.*(frequency(iq1)**2 + 
     |             frequency(io1)**2 + frequency(ik1)**2)
	aap1 = alp1/3./alq1*(1.-
     |    (frequency(ip1)-frbar1)*(frequency(iq1)-frbar1)/deno1)

	bbp1 = alp1/3./alo1*(1.-
     |    (frequency(ip1)-frbar1)*(frequency(io1)-frbar1)/deno1)

	ccp1 = alp1/3./alk1*(1.-
     |    (frequency(ip1)-frbar1)*(frequency(ik1)-frbar1)/deno1)

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
	call getenv('WAVE_PATH',WAVE_PATH)
        if (WAVE_PATH.eq.' ') WAVE_PATH='.'
	lwa=lnblnk(WAVE_PATH)
	WAVE_PATH(lwa+1:lwa+1)='/'
	lwa=lwa+1
*	print *, 'Initialisation ...please be patient ...'
*	print *,'Used environment variable :'
*	print *,'WAVE_PATH=',WAVE_PATH(1:lwa)
#ifdef ASCII_FILE 
*	print *,'Reading ASCII files'
#endif
#ifdef BIMG_FILE
*	print *,'Reading BIMG files'
#endif
#ifdef GRID_FILE
*	print *,'Reading GRID files'
#endif
	lmod=lnblnk(model)
	write (0,600) WAVE_PATH(1:lwa),model(1:lmod)
600	format ('(grenoble: loading ',a,'<',a,'>',$)
601	format (',',a,$)
602	format (')')
c
	nbext=8

c
	extended(1)='q1'               
	extended(2)='o1'               
	extended(3)='k1'               
	extended(4)='2n2' 
	extended(5)='n2'               
	extended(6)='m2'               
	extended(7)='k2'               
	extended(8)='s2' 
c              
	extended(9)='p1'               
	extended(10)='nu2'               
	extended(11)='mu2'               
	extended(12)='l2'               
	extended(13)='t2'
	extended(14)='eps2'
	extended(15)='lda2'
	extended(16)='eta2'
	extended(17)='2q1'
	extended(18)='sig1'
	extended(19)='rho1'
	extended(20)='m11'
	extended(21)='m12'
	extended(22)='xi1'
	extended(23)='pi1'
	extended(24)='phi1'
	extended(25)='tta1'
	extended(26)='j1'
	extended(27)='oo1'
 
	unit2=20
             
	nb=0      
742	continue
      do nn=1,nbext 

	  nb=nb+1 
	  lex=lnblnk(extended(nn))
	  pathname=WAVE_PATH
	  pathname(lwa+1:)=extended(nn)
	  pathname(lwa+lex+1:)='.'//model
#ifdef BIMG_FILE
	  pathname(lwa+lex+lmod+2:)='.bimg'
#endif
c
c
120	  wave(nb)=extended(nn)
	if (rep.eq.'N') goto 101 
*	print *,'trying to read ',pathname(1:lnblnk(pathname))
	write (0,601) extended(nn)(1:lex)
        call read_tide(nb,unit2,pathname,istat) 

	  if (istat.eq.0) then 
	    nbwave=nb
          do i=1,9         
            if(wave(nb).eq.name(i)) then
              freq(nb)=frequency(i)
              num(nb)=code(i)   
c              write(*,*) wave(nb),real(freq(nb)),' degree/hour'
	        freq(nb)=freq(nb)*pi/180.d+00
              goto 100
            endif
	  enddo 
	    stop 'unknown wave...'
100	    continue
	  else
	    nb=nb-1
	  endif
101	continue
	enddo  

10	format(i3,3x,a6,3x,d14.11) 
c
c nine (8) wave are read from the files 19 are inferred by admittance
c
	nbwave=27
          do i=9,nbwave   
	      nb=nb+1   
	      wave(nb)=extended(i)
              freq(nb)=frequency(i)
              num(nb)=code(i)   
c              write(*,*) wave(nb),real(freq(nb)),' degree/hour'
	        freq(nb)=freq(nb)*pi/180.d+00
	  enddo 
c	do i=1,nbwave
c	 print *,i,freq(i)/pi*180.
c	enddo
                 
	call init_corrections(0.0d0)
*	print *, 'Now ready !!!'
	write (0,602)
      return
      end
c###
      subroutine read_tide(nb,iunit,pathname,istat)
c------------------------------------------------------------------------------
c  PROGRAMME :read_tide.f
c
c  DESCRIPTION : Read data files either  ASCII formatted or BIMG files
c                (Check the CPPFLAG in Makefile)
c                to convert ASCII to BIMG use ascii2bimg.f program
c
c
c
c  PROGRAMMEURS :  F. LYARD
c 
c  DATE derniere retouche : 13/03/1995
c------------------------------------------------------------------------------
       include 'common.h'   

	integer iunit,gridrd4
	character*255 pathname,gridname
	character*80 line1,line2,line3,line4,line
                    
	real a(30),G(30)  

	spec=32767.

	istat=0
c
c default files are ASCII files, (BIMG_FILE not defined)
c
#ifdef ASCII_FILE    
	 open(iunit,file=pathname,status='old',err=100)
	read(iunit,'(a)') line
	read(line,*) xmin,xma
	line1=line(19:)
* 	print *,line1
	read(iunit,*) ymin,yma
	read(iunit,*) dx,dy
	read(iunit,*)   ni,nj
	read(iunit,*) amask,Gmask

	do j=1,nj
	  do i=1,ni,30
	    read(iunit,10) (a(k),k=1,30)
	    read(iunit,10) (G(k),k=1,30)
	    do k=1,30
	      ii=k+(i-1)
	      if (a(k).ne.amask.and.G(k).ne.Gmask) then
	        wrp(ii,j,nb)=a(k)*cos(G(k)*rad)
	        wip(ii,j,nb)=a(k)*sin(G(k)*rad)
	      else
	        wrp(ii,j,nb)=spec
	        wip(ii,j,nb)=spec
	      endif
	    enddo
	  enddo

c
c tide(-0.5,phi)=tide(359.5,phi)
c convenient for interpolation
c over a global solution
c
          wrp(0,j,nb)=wrp(ni,j,nb)
	  wip(0,j,nb)=wip(ni,j,nb)
          wrp(ni+1,j,nb)=wrp(1,j,nb)
          wip(ni+1,j,nb)=wip(1,j,nb)

	enddo	    
10	format(30f7.2)
20	format(30f7.1)

	close(iunit) 
#endif
#ifdef GRID_FILE
	l=index(pathname,' ')-1
	ni=722
	nj=361
	gridname=pathname
	gridname(l+1:)='.amp'
	ios=gridrd4(gridname,ni,nj,wrp(1,1,nb),
     |		xmin,xma,ymin,yma,zmin,zma)
	if (ios.ne.0) goto 100
	ni=722
	nj=361
	gridname(l+1:)='.pha'
	ios=gridrd4(gridname,ni,nj,wip(1,1,nb),
     |		xmin,xma,ymin,yma,zmin,zma)
	if (ios.ne.0) goto 100
	dx=(xma-xmin)/(ni-1)
	dy=(yma-ymin)/(nj-1)
	ni=ni-1
	amask=1e30
	do j=1,nj
	   wrp(0,j,nb)=wrp(ni,j,nb)
	   wip(0,j,nb)=wip(ni,j,nb)
	enddo
	 do j=1,nj
	  do i=0,ni+1
	   a(1)=wrp(i,j,nb)
	   G(1)=wip(i,j,nb)
	      if (a(1).lt.amask.and.G(1).lt.amask) then
	        wrp(i,j,nb)=a(1)*cos(G(1)*rad)
	        wip(i,j,nb)=a(1)*sin(G(1)*rad)
	      else
	        wrp(i,j,nb)=spec
	        wip(i,j,nb)=spec
	      endif
	  enddo
	 enddo	
#endif
#ifdef BIMG_FILE
	open(iunit,file=pathname,status='old',form='unformatted',
     |       err=100)
	read(iunit)line1
	read(iunit)line2
	read(iunit)line3
	read(iunit)line4
	read(iunit)nimax,njmax,nz,nt,ndim,icod
* 	print *, line1
* 	print *, line2
* 	print *, line3
* 	print *, line4
	ni=nimax-1
	nj=njmax
	read(iunit)xmin,ymin,dx,dy,amask
	read(iunit) z0
	read(iunit) tim_tag
	read(iunit) ((wrp(i,j,nb),i=0,ni),j=1,nj)
	read(iunit) ((wip(i,j,nb),i=0,ni),j=1,nj)
c row ni+1 is filled with value at i=1 for interpolation
        do j=1,nj
          wrp(ni+1,j,nb)=wrp(1,j,nb)
          wip(ni+1,j,nb)=wip(1,j,nb)
        enddo
	 do j=1,nj
	  do i=0,ni+1
	   a(1)=wrp(i,j,nb)
	   G(1)=wip(i,j,nb)
	      if (a(1).ne.amask.and.G(1).ne.amask) then
	        wrp(i,j,nb)=a(1)*cos(G(1)*rad)
	        wip(i,j,nb)=a(1)*sin(G(1)*rad)
	      else
	        wrp(i,j,nb)=spec
	        wip(i,j,nb)=spec
	      endif
	  enddo
	 enddo	
	close(iunit)
#endif
    
	ni=ni+1

*       print *,'  Wave ',wave(nb),' ok...' 
	return
            
100	continue
       print *,' ERROR opening data file', pathname
	istat=-1
	stop 'FATAL ERROR'

      return
      end
