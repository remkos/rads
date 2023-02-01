      program gravcorr
      integer lmax,j,k,l,m,nx/0/,ny/0/,lnblnk
      parameter (lmax=70)
      real*8 d,sig
      real*8 c(0:lmax,0:lmax),s(0:lmax,0:lmax)
      real*8 sigc(0:lmax,0:lmax),sigs(0:lmax,0:lmax)
      real*8 c0(0:lmax,0:lmax),s0(0:lmax,0:lmax)
      character*80 line,filenm1,filenm2

      do l=0,lmax
         do m=0,lmax
	    c(l,m)=0
	    s(l,m)=0
	    sigc(l,m)=0
	    sigs(l,m)=0
	 enddo
      enddo
	    

      call getarg(1,filenm1)
      open (10,file=filenm1,status='old')

10	 read (10,*,end=19) j,l,m,k,d,sig
	 if (l.eq.0 .and. m.eq.0) then
	 else if (k.eq.1) then
	    c(l,m)=-d
	    sigc(l,m)=-sig
	 else
	    s(l,m)=-d
	    sigs(l,m)=-sig
	 endif
	 nx=max(nx,l)
	 ny=max(ny,m)
      goto 10
19    continue
      close (10)


      call getarg(2,filenm2)
      open (10,file=filenm2,status='old')

      read (10,550) line
      write (*,550) "GRAVITY MODEL -- DGM-EXX -- ERS Tailored EGM96"
      read (10,550) line
      l=lnblnk(line)
      write (*,550) line(:l)
20    read (10,600,end=90) l,m,c0(l,m),s0(l,m)
550   format (a)
600   format (14x,2i3,9x,2d15.8)
610   format ('GCOEF 1',7x,2i3,9x,2d15.8)
      write (*,610) l,m,c0(l,m)+c(l,m),s0(l,m)+s(l,m)
      goto 20

90    close (10)

      do l=0,lmax
	 do m=0,lmax
	    if (sigc(l,m).eq.0d0) sigc(l,m)=1d35
	    if (c0(l,m).eq.0d0) c0(l,m)=1d35
	    if (sigs(l,m).eq.0d0) sigs(l,m)=1d35
	    if (s0(l,m).eq.0d0) s0(l,m)=1d35
	 enddo
      enddo
      l=index(filenm1,'.out')-1
      if (l.le.0) l=index(filenm1,' ')-1
      call gridwr8(filenm1(:l)//'_c.grd',lmax+1,lmax+1,sigc,lmax+1,
     |	0d0,dble(lmax),0d0,dble(lmax))
      call gridwr8(filenm1(:l)//'_s.grd',lmax+1,lmax+1,sigs,lmax+1,
     |	0d0,dble(lmax),0d0,dble(lmax))
*     l=lnblnk(filenm2)
*     call gridwr8(filenm2(:l)//'_c.grd',lmax+1,lmax+1,c0,lmax+1,
*    |	0d0,dble(lmax),0d0,dble(lmax))
*     call gridwr8(filenm2(:l)//'_s.grd',lmax+1,lmax+1,s0,lmax+1,
*    |	0d0,dble(lmax),0d0,dble(lmax))
      end
