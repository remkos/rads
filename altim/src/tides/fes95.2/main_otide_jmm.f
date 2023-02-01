	program main_otide_jmm
c------------------------------------------------------------------------------
c  PROGRAMME : main_otide_jmm.f
c
c  DESCRIPTION : Example of the use of the subroutine otide. This program
c               is given as a skeleton for building its own software.
c               Please NOTE that hour is real*8.
c
c
c  PROGRAMMEURS : J.M. MOLINES 
c 
c  DATE derniere retouche : 13/03/95
c  Changed to include correction for the long period tides according
c  to the programs by D. E. Cartwright, and  documented by r. ray 
c       November 1994 , Ole Baltazar Andersen, KMS 
c
c------------------------------------------------------------------------------
	real*4 lat,long
	real*8 hour,ts,lpeqtide
	integer istat
	character*20 model
c define the extension of the data files (e.g. m2.legi)
	model='fes94.1'
	model='sra95.1'
c Exemple for the use of routine otide:
c Latitudes and longitudes are to be given in decimal degrees,
c 
c Longitudes can be given either between 0 et 360 deg. or -180 et 180 deg.
c Date is given in "Jnasa" (as on AVISO/CDROM) and the time is given
c in hours, during the day (double precision, between 0.d0 and 24.d0)
c For instance Jnasa=12692 for 1 October ,1992.
c
888	continue
	read *,lat,long
	if (lat.eq.999.) goto 889
	read*,jnasa,hour
	ts=jnasa+36204+hour/24.d0
	call otide(model,tide,lat,long,jnasa,hour,istat)
	print 101,jnasa,hour,lat,long,tide*10, istat
100	format(a,f8.1)
101     format (i5,4x,f15.12,3x,f8.4,3x,f8.4, f10.2,2x,i1)
c
c The tide output is given in millimeters. It is pure ocean tide from the
c diurnal and semidiurnal band. Load tide, solid tide and long period tides
c have to be added for altimetric corrections.
c The missing value is indicated by tide=32767
c
	goto 888
889	print *,' Bye !'
	end
