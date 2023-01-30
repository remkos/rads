        program ascii2bin
c-----------------------------------------------
c
c This program convert the ASCII distributed files for
c FES Tide Model to binary files in order to save
c disk space and time in reading the files.
c Take into account compatibility with 
c LEGOS binary file (bimg format)
c
c First version : Jean-Marc Molines 15/11/1994
c Revised       : Fabien Lefevre 09/08/2000
c
c------------------------------------------------
        integer*4 nimax, njmax, iunit, i, j, k, l, m, strlen
        parameter (nimax=1441, njmax=721)
        integer*4 ni, nj, nk, nt, nd, icode, dn, ios
        real*4 xmin,ymin,xmax,ymax
        real*4 dx, dy, mask, level, pi, re, im
        real*4 wra(nimax,njmax), wrg(nimax,njmax)
        complex*8 wcomp(nimax,njmax)
        character*126 file_asc, file_bin, tabfile_asc(16)
        character*80 comment(4)

        iunit=10
        pi=acos(-1.)

        tabfile_asc(1)  = 'M2_fes99.asc'
        tabfile_asc(2)  = 'S2_fes99.asc'
        tabfile_asc(3)  = 'N2_fes99.asc'
        tabfile_asc(4)  = 'K2_fes99.asc'
        tabfile_asc(5)  = '2N2_fes99.asc'
        tabfile_asc(6)  = 'K1_fes99.asc'
        tabfile_asc(7)  = 'O1_fes99.asc'
        tabfile_asc(8)  = 'Q1_fes99.asc'
        tabfile_asc(9)  = 'M2_drfes99.asc'
        tabfile_asc(10) = 'S2_drfes99.asc'
        tabfile_asc(11) = 'N2_drfes99.asc'
        tabfile_asc(12) = 'K2_drfes99.asc'
        tabfile_asc(13) = '2N2_drfes99.asc'
        tabfile_asc(14) = 'K1_drfes99.asc'
        tabfile_asc(15) = 'O1_drfes99.asc'
        tabfile_asc(16) = 'Q1_drfes99.asc'

c--- Read input parameter
        do m=1,16
          file_asc=tabfile_asc(m)
          l=strlen(file_asc)

          print *,' Reading ASCII file ... be patient ...'
          print *,' File : ', file_asc(1:l)
          open (unit   = iunit,
     &          file   = file_asc,
     &          iostat = ios,
     &          access = 'sequential',
     &          status = 'old',
     &          form   = 'formatted',
     &          err    = 100)

c--- Read ascii file
          read(iunit,*,err=101) xmin,xmax
          read(iunit,*,err=101) ymin,ymax
          read(iunit,*,err=101) dx,dy
          read(iunit,*,err=101) ni,nj
          read(iunit,*,err=101) mask,mask

c          write(6,'(2f9.3)') xmin,xmax
c          write(6,'(2f9.3)') ymin,ymax
c          write(6,'(2f9.3)') dx,dy
c          write(6,'(2i6)')   ni,nj
c          write(6,'(2f9.3)') mask

          do j=1,nj
            do k=1,ni,30
              read(iunit,*) (wra(i,j),i=k,MIN(ni,k+29))
              read(iunit,*) (wrg(i,j),i=k,MIN(ni,k+29))
            enddo
          enddo

          close(iunit) 

c--- Write binary file : compatibility with Bimg LEGOS format
          l=strlen(file_asc)
          file_bin=file_asc(1:l-4)//'.bin'

          print *,' Writing binary file ... will not be so long ...'
          open (unit = iunit,
     &          file=file_bin,
     &          form='unformatted',
     &          access='sequential',
     &          status= 'unknown',
     &          err    = 102)

          comment(1)='Finite Element Solution 99 (FES99)'
          comment(2)='Copyright : LEGOS/GRGS 1999'
          comment(3)='14, av Edouard Belin'
          comment(4)='31400 Toulouse - FRANCE'
          do i=1,4
            write(iunit, err=102) comment(i)
          enddo

          nk=1
          nt=1
          nd=1
          icode=1
          write(iunit) ni,nj,nk,nt,nd

c--- Shift 0 to 360 degrees
          write(iunit, err=102) 0,ymin,dx,dy,mask,icode
          level=1.
          write(iunit, err=102) level
          write(iunit, err=102) level

c--- Shift 0 to 360 degrees
          dn=-xmin/dx
          do i=1,ni
            if(i.le.dn) then
              k=i+dn+1
            else
              k=i-dn
            endif
            do j=1,nj
              if((wra(i,j).eq.mask).or.(wrg(i,j).eq.mask)) then
                wcomp(k,j)=cmplx(mask,mask)
              else
                re=wra(i,j)*cos(-wrg(i,j)*pi/180.)/100.
                im=wra(i,j)*sin(-wrg(i,j)*pi/180.)/100.
                wcomp(k,j)=cmplx(re,im)
              endif
            enddo
          enddo
          write(iunit, err=102) ((wcomp(i,j),i=1,ni),j=1,nj)

          close(iunit)
          
          goto 110
 100      print *,' Problem opening ASCII file ... please check'
          goto 110
 101      print *,' Problem opening header file ... please check'
          goto 110
 102      print *,' Problem writing BINARY file ... please check'
 110      continue
 
        enddo
        end
        
        
 
