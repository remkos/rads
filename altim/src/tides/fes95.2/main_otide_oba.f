	program main_otide_oba
c------------------------------------------------------------------------------
c  PROGRAMME : main_otide_oba.f
c
c  DESCRIPTION : Example of the use of the subroutine otide. This program
c               is given as a skeleton for building its own software.
c               Please NOTE that hour is real*8.
c
c
c  PROGRAMMEURS : J.M. MOLINES, O.B. ANDERSEN
c
c  DATE last touch : 13/03/95
c
c  Changed to include correction for the long period tides according
c  to the programs by D. E. Cartwright, and  docomented by r. ray 
c       November 1994 , Ole Baltazar Andersen, KMS 
c
c------------------------------------------------------------------------------
	real*4 lat,long
	real*8 heure,dt,mjul85,ts,tlp,allt
        character*70 fname,oname
	character*20 model
	
c
c
	print *, 'Tide correction program from GRENOBLE.'
	print *, 'The program adds ocean tide correction or'
	print *, 'elastic ocean tide correstion to altimetry'
	print *, 'using Grenoble model or Andersen 1995.1 adjusted model'
	print *, ' '
	print *, 'NOTE: Long period correction included using equilibrium'
	print *, '      response as implemented by Carwright & Ray'
	print *, '      For the Andersen 1995.1 model a Mediterranean '
	print *, '      ocean tide correction is included using a '
	print *, '      0.5 degree re-interpolated Canceill model '
	print *, ' '
	print *, '                          Ole B. Andersen 19/12/1994'
	print *, ' '
	print *,' Another example is given using different time format'
	print *,' in the program main_otide_jmm.f'
	print *,'                           Jean-Marc Molines 13/03/1995'
c
c The variable model contains the extension of the data files. For instance
c if your data files are m2.fes95.1, s2.fes95.1 etc, model='fes95.1' 
c If you decided to use binary files, the files wille be m2.fes95.1.bimg, BUT
c model is still 'fes95.1'
c
	model='fes95.1'
c
c.....INPUTS:
c       LAT  - Latitude(degrees)
c       LON  - Longitude(degrees)
c       TIME - Time after  1: Jan 1, 0 hrs, 1992(seconds)
c       TIME - Time after  2: Jan 1, 0 hrs, 1985(seconds)
c.....OUTPUTS:
c       TIDE - Ocean tide correction in meters
c              TIDE = 0 indicates that no ocean tide correction is
c              applied (no coverage of model)

        print *,'Enter filename for input data file '
        print *,'Default data format: time,lat,long,.......'

        read(*,'(A)') fname
        print *,'Enter output filename '
        read(*,'(A)') oname
        print *,'Please enter reference time frame (1 or 2)'
        print *,'1: Ref time 1/1 - 1992,0 hour (Default TOPEX/POSEIDON)' 
        print *,'2: Ref time 1/1 - 1985,0 hour (Default ERS1 or GEOSAT)' 
        read(*,*) itime
        if (itime.lt.1.or.itime.gt.2) stop 'Value must be equal 1 or 2'

        open(10,file=fname,err=9003)
c unit 20 is opened and used in the package
        open(25,file=oname,err=9003)
         
 9001   continue
        read (10,*,end=9002) dt,lat,long
        if (itime.eq.1) dt = dt+2556*86400.0
        jday = idint(dt/86400.0)
        jnasa = 9862+jday
        heure = (dt-86400.0*jday)/3600.0

c       write(21,*) jday,jnasa,heure
c       if (lat.eq.999.) goto 889
c       print *,' Check for correct time for TOPEX : Jnasa , heure '
c       print *,' This should correcspont to the TOPEX record:'
c       print *,' Jnasa = Tim_Moy_1'
c       print *,' heure = (Tim_Moy_2/1000.d0+Tim_Moy_3/1.d6)/3600.d0'
c       read*,jnasa,heure

	call otide(model,tide,lat,long,jnasa,heure,istat)
	print *,tide,lat,long,jnasa,heure,istat
	tide=tide/100.

c------------------------------------------------------------------- 
c istat = number of points used in the interpolation
c istat = 0 the point (rlon,rlat) is out of domain
c istat = 1,2,3 the point (rlon,rlat) is on the border of the domain
c oba: might be adequate for tidal interpolation towards the coast
c istat = 4 the point (rlon,rlat) is fully inside the domain
c-------------------------------------------------------------------
      if (abs(tide).gt.9.0) tide = 999.999
      if (istat.lt.4)       tide = 999.999

C Correction for low frequency tides - using R. Rays approach.

C MODIFIED JULIAN DAY AT START 1985 IN SECONDS     
      MJUL85 =  46066.D0*24.D0*3600.D0
C
C CONVERT TO MODIFIED JULIAN DATE. 
c MJD OF 0 HOURS 1 JANUARY 1985 IS 46066.
      TS=DT+MJUL85           
C                                                                     
C NOTE THAT LPEQMT IS ONLY AN APPROXIMATE (BUT FAST) ESTIMATE OF THE 
C LONG PERIOD OCEAN TIDE BUT SHOULD BE OK FOR ALTIMETRY  
C                      
      CALL LPEQMT( TS, LAT, TLP )        

C ADD ON THE SMALL LONG PERIOD BIT TO GET THE 'TOTAL' TIDE 

      ALLT=TIDE+TLP/100.0
      if (allt.gt.990.0) allt = 999.999

c     Invoke the following line if output is wanted in centimeters
c     ALLT = ALLT*100.0

      write(25,557) DT,LAT,LONG,ALLT
 557    format(f12.1,2f10.4,f10.3)
	
        goto 9001
 9002   continue
        goto 9004
 9003   continue
         write(*,*) ' Could not find input/output file names '
 9004   continue

889	print *,' Bye !'
        close(10)
        close(25)
	end

