      SUBROUTINE PRERROR (NMEAS)

* This subroutine prints a overview of the various errors concerning
* the translated data.

      integer i,j,nmeas(10,*),ntot(10)/10*0/

600   format (//'Overview of measurement errors:'/
     |' A = Epoch time scale unknown (measurement deleted)'/
     |' B = BIPM time correction required before start of times.data',
     |' table (1st value used)'/
     |' C = BIPM time correction required after end of times.data',
     |' table (last value used)'/
     |' D = GPS time correction required before start of times.data',
     |' table (zero value used)'/
     |' E = GPS time correction required after end of times.data',
     |' table (last value used)'/
     |' F = Measurement not corrected for CG offset (unknown',
     |' satellite)'/
     |' G = No proper wavelength could be found in system.data',
     |' table (600 nm assumed)'/
     |' H = Format error (measurement deleted)'/
     |' I = Checksum error (ignored)'//
     |'Station Total     A     B     C     D     E     F     G     H',
     |'     I')
      do i=1,9999
	 do j=2,10
	    nmeas(1,i)=nmeas(1,i)+nmeas(j,i)
	    ntot(j)=ntot(j)+nmeas(j,i)
	 enddo
	 ntot(1)=ntot(1)+nmeas(1,i)
      enddo
      if (ntot(1).gt.0) then
         write (*,600)
         do i=1,9999
            if (nmeas(1,i).gt.0)
     |		write (*,'(i7,10i6)') i,(nmeas(j,i),j=1,10)
         enddo
         write (*,'("  Total",10i6)') ntot
      endif
      end
