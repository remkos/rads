c
c     this program selects a subset of the JGM-2 normal
c     equations, and writes in E-matrix format
c
c     UNIT 10: CRS MATRIX
c     UNIT 20: SELECTED COEFFICIENTS
c     UNIT 30: E-MATRIX
c
      implicit real*8(a-h,o-z)
      parameter(nparam=5035)
      parameter(maxcoef=2000)
      parameter(maxdeg=70)
 
      character*7 char
      
      real*8 norm(maxcoef*(maxcoef+1)/2)
      real*8 header(30),label(3,nparam)
      real*8 solveh(nparam),sigmah(nparam),scaleh(nparam)
      real*8 cnm(0:maxdeg,0:maxdeg,2)
      
      real row(nparam)
      real solve(nparam),sigma(nparam),scale(nparam)

      integer ics(nparam),iord(nparam),ideg(nparam)
      integer icsh(maxcoef),iordh(maxcoef),idegh(maxcoef)
      integer itest(0:maxdeg,0:maxdeg,2)
      integer index(0:maxdeg,0:maxdeg,2)
      
      data header/-931234567899.d0,-931234567899.d0,979.0,
     c  0.,100.0,100.0,100.0,1.0,1.0,300.0,0.,0.,0.,9208.08d0,
     c  1.0,299792458.d0,14*0./

      open(unit=10,status='old',form='unformatted',file=
     c    '/u2a/pieter/Covar/JGM3.70.NOR')
      open(unit=20,status='old',form='formatted',file=
     c   'gcoef.out')
c     open(unit=30,status='new',form='formatted',file=
      open(unit=30,status='unknown',form='unformatted',file=
     c   '/u2a/geodyn/gdyn9404/ers2/emat/solve/JGM3.CONSTR')
      
      read(10) n
      do 10 i=1,n
      read(10) ics(i),ideg(i),iord(i),solve(i),sigma(i),scale(i)
      solve(i)=solve(i)/scale(i)
   10 continue
      
      do 20 l=0,maxdeg
      do 20 m=0,l
      do 20 k=1,2
      itest(l,m,k)=0
   20 continue
      
      number=0
   30 read(20,300,end=40) char,l,m,coefc
  300 format(a7,7x,2i3,4x,d20.8,15x,d13.1)
      ic=1
      if (char.eq.'GCOEFS1') ic=2
      itest(l,m,ic)=1
      number=number+1
      index(l,m,ic)=number
      icsh(number)=1
      iordh(number)=m
      idegh(number)=l
      cnm(l,m,ic)=coefc
c     
c     determine label
      if (ic.eq.1) then
      label(1,number)=5030010.d7+l*1.d4+m*1.d0
      else
      label(1,number)=5030020.d7+l*1.d4+m*1.d0
      endif
      goto 30
      
   40 continue
      write(*,*) number
      
      do 50 i=1,n
      read(10) (row(j),j=1,i)
      if (itest(ideg(i),iord(i),ics(i)).eq.1) then
      index1=index(ideg(i),iord(i),ics(i))
      solveh(index1)=cnm(ideg(i),iord(i),ics(i))
      sigmah(index1)=sigma(i)
      scaleh(index1)=scale(i)
      do 60 j=1,i  
      if (itest(ideg(j),iord(j),ics(j)).eq.1) then
      index2=index(ideg(j),iord(j),ics(j))
      if (index1.gt.index2) then
      norm(index1*(index1-1)/2+index2)=row(j)/scale(i)/scale(j)
      else
      norm(index2*(index2-1)/2+index1)=row(j)/scale(i)/scale(j)
      endif
      endif
   60 continue
      endif
   50 continue
      
c
c     adjust header record
      header(3)=number
      write(30) header
c     write(30,'(6d24.17)') header
      
c
c     write parameter group identifier
      pgroup=5.d13
      write(30) (pgroup,i=1,nint(header(15)))
c     write(30,'(6d24.17)') (pgroup,i=1,nint(header(15)))

c
c     write number of parameters record
      dnum=number*1.d0
      write(30) (dnum,i=1,nint(header(15)))
c     write(30,'(6d24.17)') (dnum,i=1,nint(header(15)))
      
c
c     write labels
      do 80 i=1,number
      label(2,i)=0.d0
      label(3,i)=0.d0
   80 continue
      dword=dfloat(3*number+1)
      dresid=2.d14
      write(30) dword,(label(1,j),j=1,number),
     c (label(2,j),j=1,number),(label(3,j),j=1,number),dresid
c     write(30,'(6d24.17)') dword,(label(1,j),j=1,number),
c    c (label(2,j),j=1,number),(label(3,j),j=1,number),dresid

c     
c     write parameters values record (=a priori)
      write(30) dnum,(solveh(i),i=1,number)
c     write(30,'(6d24.17)') dnum,(solveh(i),i=1,number)
    
c
c     write parameter a priori values record
      write(30) dnum,(solveh(i),i=1,number)
c     write(30,'(6d24.17)') dnum,(solveh(i),i=1,number)

c
c     write normal equations
      do 90 i=1,number
      drow=i*1.d0
      delem=i*1.d0
      resid=0.d0
      write(30) drow,delem,(norm(j*(j-1)/2+i),
     c  j=i,number),resid
c     write(30,'(6d24.17)') drow,delem,(norm(j*(j-1)/2+i),
c    c  j=i,number),resid
   90 continue
 
      end
