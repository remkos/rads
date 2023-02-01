      PROGRAM dorisres
* 
* Processes GEODYN residual files for DORIS data analysis and editing
*
*-----------------------------------------------------------------------
      implicit none 

* Declaration of variables

      integer maxobs,maxpas
      parameter(maxobs=200,maxpas=2500)

      real*8 times(maxobs,maxpas),resid(maxobs,maxpas),
     |sigma(maxobs,maxpas),deriv(maxobs,maxpas),elev(maxobs,maxpas)
      integer ipass(4,maxpas)
      character cpass(maxpas)*8,filenm*80
      logical timesort
      integer*4 i, j, iarg

* Scan the arguments

10    iarg=iarg+1
      call getarg(iarg,filenm)
      if (filenm.eq.'-t') then
	 timesort=.true.
	 goto 10
      else if (filenm.eq.' ') then
         filenm='fort.19'
      endif

* Open residual file and read it

      call resread(filenm,maxobs,maxpas,ipass,cpass,
     |	times,resid,sigma,deriv,elev)

      do i=1, maxpas
        if(ipass(4,i).eq.40) then
          print *, ''
          print *, ''
          print *, '#', cpass(i)
          do j=1, ipass(3,i)
            print *, times(j,i),resid(j,i),elev(j,i)
          enddo
        endif
      enddo

      end


***********************************************************************
* Subroutine meanrms: Compute mean and RMS of a number of data values

      subroutine meanrms(n,values,mean,rms)
      integer n, i
      real*8 mean,rms
      real*8 values(*)
      if (n.gt.0) then
       do i=1,n
         mean = mean + values(i)
         rms = rms + values(i)**2
       enddo
      endif
      mean=mean/n
      rms=sqrt(rms/n)
      end
