**GLOBPRES - Determine global mean pressure over oceans
*+
      SUBROUTINE GLOBPRES (TYPE, UTC, VALUE)
      INTEGER*4 TYPE
      REAL*8    UTC, VALUE

* Determine the global mean pressure over oceans by linear
* interpolation in a table of 6-hourly global mean pressure
* values based on ECMWF and NCEP grids of sea surface pressure.
*
* The value returned can be the one based on NCEP grids
* (available for the GEOSAT and GFO time frames) or on
* ECMWF grids (available for the T/P time frame). Also
* smoothing can be applied. The smoother is a windowed sinc
* function with a (2 days)**(-1) cut-off frequency.
*
* Input arguments:
*  TYPE  : Type of value requested
*      =1: Value from ECMWF pressure fields
*      =2: Idem, smoothed
*      =3: Value from NCEP pressure fields
*      =4: Idem, smoothed
*  UTC   : Time in UTC seconds since 1.0 Jan 1985
*
* Output argument:
*  VALUE : Global mean pressure over oceans (mbar)
*-
* 26-Feb-2002 -- Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer unit,freeunit,nval/-1/,i,nmax
      parameter (nmax=100000)
      real*8  t,t0/1d30/,t1/-1d30/,press(4,nmax),dt,x
      character*80 filename,line
      save t0,dt,press,nval

* On first call: open table with pressure values. Load them in memory.

      if (nval.lt.0) then
         unit=freeunit()
	 filename='/user/altim'
         call checkenv('ALTIM',filename,i)
         filename(i+1:)='/data/tables/global_pressure.txt'
         nval=0
         open (unit,file=filename)
10       read (unit,'(a)',end=100) line
         if (line(:1).eq.'#') goto 10
         nval=nval+1
	 if (nval.gt.nmax) call fin('GLOBPRES: too many values')
         read (line,*) t,(press(i,nval),i=1,4)
	 t0=min(t0,t)
	 t1=max(t1,t)
         goto 10
100      close (unit)
	 t0=(t0-46066d0)*86400d0
	 t1=(t1-46066d0)*86400d0
         dt=(t1-t0)/(nval-1)
      endif

      x=(utc-t0)/dt+1
      i=int(x)
      if (i.lt.1 .or. i.ge.nval) then
         value=1013.3d0
      else
         x=x-i
         value=press(type,i)*(1-x)+press(type,i+1)*x
      endif
      end
