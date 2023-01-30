**BLTWEEN -- Quick search for measurement with timelabel TIMX
*+
      subroutine bltween(iblk,rdata,adata,
     |                  meas0,measlo,meashi,timelo,timehi,timx)
************************************************************************
* Seek meas0 : last  measurement where t<timx                          *
*      meas1 : first measurement where t>=timx                         *
* Modification : Does not jump to the midpoint of the interval but     *
*                computes roughly where x-over can be in file          *
************************************************************************
      implicit none
      integer*4 iblk,meas,meas0,meas1,measlo,meashi
      integer*2 adata(3,2)
      real*8 rdata(5,2),tobs,timx,dt,timelo,timehi,dtmin
      parameter (dtmin=979921d-6)

* Make a first guess of the range of measurement numbers (meas0,meas1).
* Check for data beyond the track limits.

      dt=timehi-timx
      if (dt.lt.dtmin) then
	 meas0=meashi-1
         meas1=meashi
	 goto 100
      else
	 meas0=max(measlo,meashi-int(dt/dtmin)-1)
      endif
      dt=timx-timelo
      if (dt.lt.dtmin) then
	 meas0=measlo
         meas1=measlo+1
	 goto 100
      else
         meas1=min(meashi,measlo+int(dt/dtmin)+1)
      endif

* Guess the observation closest to the crossover to be half way between
* the beginning and the end of the current range.

10    continue
      meas=(meas0+meas1)/2
      call blread(iblk,meas,rdata(1,1),adata(1,1))
      tobs=rdata(1,1)

* If the found crossover is before the crossover we seek then the
* crossover we seek is between the found crossover (index O) and the
* previously found upper bound (index 1). When the time difference
* between the found crossover and the crossover we seek is smaller than
* a threshold level we can safely assume we've found it. The found
* crossover cannot be outside the 0-1 range so we don't need to check for
* that.
      if (tobs.lt.timx) then
         meas0=meas
         if ((timx-tobs).lt.dtmin) meas1=meas+1
      else
         meas1=meas
         if ((tobs-timx).lt.dtmin) meas0=meas-1
      endif

      if (meas1-meas0.gt.1) goto 10

100   continue
      call blread(iblk,meas0,rdata(1,1),adata(1,1))
      call blread(iblk,meas1,rdata(1,2),adata(1,2))
      end
