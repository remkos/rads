**GIIESUM2 -- Process Geodyn IIE summary files (giie.sum)
*+
      program giiesum2

* This program processes Geodyn IIE (GIIE) summary files (giie.sum) created
* from the GIIE standard output (giie.out) by the utility giiesumm.
*
* Features of the current version:
* - giie.sum content is read from standard input. This may be a
*   concatination of giie.sum files.
* - Solved-for coordinates are read and stored per statio1. All
*   coordinates are converted to the latest coordinate epoch using the
*   latest station velocities.
* - A combined coordinate solution is created for each station based on
*   the standard deviations and correlations of XYZ coordinates.
* - One file per station is created (ssss.sum, where ssss is the station
*   ID) containing: all converted station coordinates plus error
*   information, the combined solution plus error information.
*
*-
* 24-Jan-2000 - Created by Remko Scharroo
* 31-Jan-2000 - Additional output in genstapos format
* 14-Mar-2001 - Changed addition of constraint
*-----------------------------------------------------------------------
      implicit none
      integer ista,i,mrec,nrec,irec
      parameter (mrec=10000)
      character line*160,text*32
      integer sta(mrec),h(3),kp(3),mdate,nsolve,nr(1000:9999)
      integer statinfo,lnblnk
      real*8 pos(20,mrec),vel(4,1000:9999),xdate,small,edit
      real*8 ecc(3),dpos(3),r,lat,lon,height
      real*8 mean(3,1000:9999),rms(3,1000:9999),ref(6,1000:9999)
      real*8 ata(6),ata0(6),ata1(6),atb(3),atb1(3),v(3),d,epoch,sec85
      parameter (small=9d-13,edit=3d0)
      logical used(mrec)

* Initialize

      call matsy1(3,h)
      nrec=0
      do ista=1000,9999
	 nr(ista)=0
	 do i=1,3
	    mean(i,ista)=0d0
	    rms(i,ista)=0d0
	 enddo
      enddo
      do i=1,6
         ata0(i)=0
      enddo

* Store station position and error matrices
* POS(i,nrec) contains:
* - I= 1, 3 : Adjusted coordinates (X,Y,Z)
* - I= 4, 6 : Coordinate error estimates (sigma X, sigma Y, sigma Z)
* - I= 7, 9 : Error covariances (rho X-Y, rho X-Z, rho Y-Z)
* - I=10,12 : Apriori coordinates (X', Y', Z')
* - I=13,15 : Apriori error estimates (sigma' X, sigma' Y, sigma' Z)
* - I=16,18 : Error covariances (rho' X-Y, rho' X-Z, rho' Y-Z)
* - I=19    : Coordinate epoch
* - I=20    : Solution epoch
* VEL(i,ista) contains for the LAST occurance of station ista:
* - I= 1, 3 : Velocity components (vel X, vel Y, vel Z)
* - I= 4    : Coordinate epoch

10    read (*,550,end=200) line
      if (line(:6).eq.'STA   ') then
         nrec=nrec+1
	 if (nrec.gt.mrec) call fin('Too many records')
         read (line(7:),*) ista,(pos(i,nrec),i=1,9)
	 sta(nrec)=ista
	 pos(20,nrec)=d
      else if (line(:6).eq.'STAAPR') then
	 read (line(7:),*) ista,(pos(i,nrec),i=10,18)
      else if (line(:6).eq.'STAVEL') then
	 read (line(7:),*) ista,(vel(i,ista),i=1,4)
	 vel(4,ista)=mdate(2,nint(vel(4,ista)))
	 pos(19,nrec)=vel(4,ista)
	 epoch=vel(4,ista)
      else if (line(:6).eq.'EPOCH ') then
         read (line(7:),*) d
	 d=sec85(4,d)/86400d0+46066
      endif
      goto 10
200   continue

* Correct station coordinates for station velocity and eccentricity
* Save reference station coordinates and sigmas

      do irec=1,nrec
	 ista=sta(irec)
	 i=statinfo(pos(20,irec),ista,text,line,line,line,line,line,ecc)
	 if (i.ne.0) then
	    write (*,*) ista,pos(20,irec)
	    call fin('Station not found in system.data')
	 endif
	 call xyzgeo(pos(1,irec),r,lat,lon,height)
	 dpos(1)=-ecc(1)*sin(lat)*cos(lon)-ecc(2)*sin(lon)
     |	         +ecc(3)*cos(lat)*cos(lon)
	 dpos(2)=-ecc(1)*sin(lat)*sin(lon)+ecc(2)*cos(lon)
     |	         +ecc(3)*cos(lat)*sin(lon)
	 dpos(3)=+ecc(1)*cos(lat)+ecc(3)*sin(lat)
	 d=(epoch-pos(19,irec))/365.25d0
         do i=1,3
	    pos(  i,irec)=pos(   i,irec)-dpos(i)+vel(i,ista)*d	! position
	    pos(9+i,irec)=pos( 9+i,irec)-dpos(i)+vel(i,ista)*d	! reference position
	    ref(  i,ista)=pos( 9+i,irec)	! save reference position
	    ref(3+i,ista)=pos(12+i,irec)	! save reference position sigma
	 enddo
      enddo

* Reduce coordinates to reference coordinates to improve computational
* efficiency.

      do irec=1,nrec
         ista=sta(irec)
	 do i=1,3
	    pos(  i,irec)=pos(  i,irec)-ref(i,ista)
	    pos(9+i,irec)=pos(9+i,irec)-ref(i,ista)
	 enddo
      enddo

* Make quick path through data. Compute mean and rms-about-mean for editing.
* Do not take into account those far away from the reference coordinates
* or those that have a sigma equal to the apriori.

      do irec=1,nrec
	 ista=sta(irec)
         d=0
	 do i=1,3
	    d=d+(pos(i,irec)/pos(12+i,irec)/edit)**2
         enddo
	 if (d.le.1 .and. pos(4,irec).lt.pos(13,irec)) then
	    do i=1,3
	       mean(i,ista)=mean(i,ista)+pos(i,irec)
	       rms(i,ista)=rms(i,ista)+pos(i,irec)**2
	    enddo
	    nr(ista)=nr(ista)+1
	 endif
      enddo
      do ista=1000,9999
         if (nr(ista).gt.1) then
	    do i=1,3
	       mean(i,ista)=mean(i,ista)/nr(ista)
	       rms(i,ista)=sqrt(rms(i,ista)/nr(ista)-mean(i,ista)**2)
	       write (*,*) ista,i,nr(ista),mean(i,ista),rms(i,ista)
	    enddo
	 endif
      enddo

* Edit records

      do ista=1000,9999
         nr(ista)=0
      enddo
      do irec=1,nrec
	 ista=sta(irec)
         d=0
	 do i=1,3
	    d=d+((pos(i,irec)-mean(i,ista))/rms(i,ista)/edit)**2
         enddo
	 if (d.le.1 .and. pos(4,irec).lt.pos(13,irec)) then
	    used(irec)=.true.
	    nr(ista)=nr(ista)+1
	 else
	    used(irec)=.false.
	 endif
      enddo

* Open file for all coordinates (input to genstapos)

      open (20,file='all.sum')
      write (20,601) xdate(epoch)

* Process station by station

      do ista=1000,9999
	 if (nr(ista).gt.1) then

* Initialize matrices

	    do i=1,6
	       ata(i)=0
	    enddo
	    do i=1,3
	       atb(i)=0
	    enddo
	    nsolve=0

* Open output file for this station


* Process all relevant records

	    do irec=1,nrec
	       if (sta(irec).eq.ista .and. used(irec)) then
		  nsolve=nsolve+1

* Build up ATAinv matrix for this record

		  ata1(1)=pos(4,irec)**2			! XX
		  ata1(3)=pos(5,irec)**2			! YY
		  ata1(6)=pos(6,irec)**2			! ZZ
		  ata1(2)=pos(7,irec)*pos(4,irec)*pos(5,irec)	! XY
		  ata1(4)=pos(8,irec)*pos(4,irec)*pos(6,irec)	! XZ
		  ata1(5)=pos(9,irec)*pos(5,irec)*pos(6,irec)	! YZ

* Do the same for the apriori ATAinv matrix for this record

		  ata0(1)=pos(13,irec)**2			! XX
		  ata0(3)=pos(14,irec)**2			! YY
		  ata0(6)=pos(15,irec)**2			! ZZ
		  ata0(2)=pos(16,irec)*pos(13,irec)*pos(14,irec)	! XY
		  ata0(4)=pos(17,irec)*pos(13,irec)*pos(15,irec)	! XZ
		  ata0(5)=pos(18,irec)*pos(14,irec)*pos(15,irec)	! YZ

* Invert ATAinv matrix for this record (creating ATA)
* Add ATA matrix for this record to the global ATA matrix
* WARNING: When inversion is imposible (if<>0 or d<small) skip the
*          processing of this record

		  call invcho(ata1,h,v,kp,3,d,i)
		  if (i.ne.0 .or. d.lt.small) then
		     write (*,*) 'ata1',d,i
		     nsolve=nsolve-1
		     goto 300
		  endif
		  call invcho(ata0,h,v,kp,3,d,i)
		  if (i.ne.0 .or. d.lt.small) then
		     write (*,*) 'ata0',d,i
		     nsolve=nsolve-1
		     goto 300
		  endif
		  do i=1,6
		     ata(i)=ata(i)+ata1(i)
		  enddo

* Multiply ATA matrix with coordinate solutions (creating ATB)
* Update global ATB matrix too.

		  call matsmv(3,h,ata1,pos(1,irec),atb1)
		  do i=1,3
		     atb(i)=atb(i)+atb1(i)
		  enddo

* Remove constraint from global ATA and ATB matrices

		  do i=1,6
		     ata(i)=ata(i)-ata0(i)
		  enddo
		  call matsmv(3,h,ata0,pos(10,irec),atb1)
		  do i=1,3
		     atb(i)=atb(i)-atb1(i)
		  enddo

* Restore reference coordinates

		  do i=1,3
		     pos(  i,irec)=pos(  i,irec)+ref(i,ista)
		     pos(9+i,irec)=pos(9+i,irec)+ref(i,ista)
		  enddo

* Print out this record (open output file)
* Print apriori coordinates for first record only

		  if (nsolve.eq.1) then
		     write (*,610) 'STAAPR',
     |			ista,(pos(i,irec),i=10,18),xdate(epoch)
		     write (line,'(i4.4,".sum")') ista
		     open (10,file=line)
		     i=statinfo(pos(20,irec),ista,text,line,line,line,line,line,ecc)
		     i=lnblnk(text)
		     write (10,600) text(:i),ista
		     write (10,620) (pos(i,irec),i=10,18),xdate(epoch),ista
		     write (10,550)
		     write (10,550)
		     write (10,550) '# Adjusted coordinates'
		  endif
		  write (10,620) (pos(i,irec),i=1,9),xdate(pos(20,irec)),ista
300		  write (*,610) 'STA   ',
     |		ista,(pos(i,irec),i=1,9),xdate(pos(20,irec))

	       endif
	    enddo 

* Add contraint to ATA matrix (diagonal elements only)

	    do i=1,3
	       ata(i+h(i))=ata(i+h(i))+ref(3+i,ista)**(-2)
	    enddo

* Invert global ATA matrix and solve global coordinate solution

	    if (nsolve.ge.1) then
	       call invcho(ata,h,v,kp,3,d,i)
	       call matsmv(3,h,ata,atb,atb1)

* Compute sigmas and correlations

	       ata1(1)=sqrt(ata(1))
	       ata1(2)=sqrt(ata(3))
	       ata1(3)=sqrt(ata(6))
	       ata1(4)=ata(2)/ata1(1)/ata1(2)
	       ata1(5)=ata(4)/ata1(1)/ata1(3)
	       ata1(6)=ata(5)/ata1(2)/ata1(3)

* Restore reference coordinates

	       do i=1,3
	          atb1(i)=atb1(i)+ref(i,ista)
	       enddo

* Print out solution

	       write (10,550)
	       write (10,550)
	       write (10,550) '# Final coordinate solution'
	       write (10,620) atb1,ata1,xdate(epoch),ista
	       write (*,610) 'SOLVED',ista,atb1,ata1,xdate(epoch)
	       close (10)
	       write (20,630) text(:6),ista,atb1
	    endif
         endif
      enddo

550   format(a)
600   format ('# Coordinate solutions: ',a,' (',i4.4,')'/
     |'# - Coordinate X, Y, Z (m)'/'# - Sigma X, Y, Z (m)'/
     |'# - Correlation XY, XZ, YZ (-)'/'# - Station ID'/'#'/
     |'# Apriori coordinates')
601   format ('# This file was created by GIIESUM2'/
     |'EPOCH     ',f15.4)
610   format(a6,i12,1x,3f18.6,7f12.3)
620   format(3f18.6,7f12.3,i6)
630   format(a6,i4.4,3f15.4/'NUVEL')
      end

      function xdate(x)
      real*8 x,xdate
      integer i,mdate
      i=int(x)
      xdate=mdate(1,i)+(x-i)
      end
