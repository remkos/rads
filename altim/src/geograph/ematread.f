      implicit real*8(a-h,o-z)
      parameter(nparam=3715)
      parameter(maxdeg=70)
      character*80 arg
      
      real*8 norm(nparam*(nparam+1)/2)
      real*8 header(30),label(3,nparam)
      real*8 solveh(nparam),sigmah(nparam),scaleh(nparam)
      
*     real row(nparam)
*     real solve(nparam),sigma(nparam),scale(nparam)

*     integer ics(nparam),iord(nparam),ideg(nparam)
*     integer icsh(nparam),iordh(nparam),idegh(nparam)
*     integer itest(0:maxdeg,0:maxdeg,2)
*     integer index(0:maxdeg,0:maxdeg,2)
      integer pgroup
      
      data header/-931234567899.d0,-931234567899.d0,979.0,
     |  0.,100.0,100.0,100.0,1.0,1.0,300.0,0.,0.,0.,9208.08d0,
     |  1.0,299792458.d0,14*0./

      call getarg(1,arg)
      open (20,file=arg,form='unformatted')

*     read(10) n
*     do 10 i=1,n
*     read(10) ics(i),ideg(i),iord(i),solve(i),sigma(i),scale(i)
*     solve(i)=solve(i)/scale(i)
*  10 continue
*     
*     do 20 l=0,maxdeg
*     do 20 m=0,l
*     do 20 k=1,2
*     itest(l,m,k)=0
*  20 continue
*     
*     number=0
*  30 read(20,300,end=40) char,l,m,coefc
* 300 format(a7,7x,2i3,4x,d20.8,15x,d13.1)
*     ic=1
*     if (char.eq.'GCOEFS1') ic=2
*     itest(l,m,ic)=1
*     number=number+1
*     index(l,m,ic)=number
*     icsh(number)=1
*     iordh(number)=m
*     idegh(number)=l
*     cnm(l,m,ic)=coefc
*     
*     determine label
*     if (ic.eq.1) then
*     label(1,number)=5030010.d7+l*1.d4+m*1.d0
*     else
*     label(1,number)=5030020.d7+l*1.d4+m*1.d0
*     endif
*     goto 30
*     
*  40 continue
*     write(*,*) number
*     
*     do 50 i=1,n
*     read(10) (row(j),j=1,i)
*     if (itest(ideg(i),iord(i),ics(i)).eq.1) then
*     index1=index(ideg(i),iord(i),ics(i))
*     solveh(index1)=cnm(ideg(i),iord(i),ics(i))
*     sigmah(index1)=sigma(i)
*     scaleh(index1)=scale(i)
*     do 60 j=1,i  
*     if (itest(ideg(j),iord(j),ics(j)).eq.1) then
*     index2=index(ideg(j),iord(j),ics(j))
*     if (index1.gt.index2) then
*     norm(index1*(index1-1)/2+index2)=row(j)/scale(i)/scale(j)
*     else
*     norm(index2*(index2-1)/2+index1)=row(j)/scale(i)/scale(j)
*     endif
*     endif
*  60 continue
*     endif
*  50 continue
*     
*
*     adjust header record
      read(20) header
      write(*,'(6d24.17)') header
      
*
*     write parameter group identifier
      read(20) pgroup
      write(*,'(6d24.17)') pgroup

*
*     write number of parameters record
      read(20) dnum
      write(*,'(6d24.17)') dnum
      number=nint(dnum)
      
*
*     write labels
      read(20) dword,(label(1,j),j=1,number),
     | (label(2,j),j=1,number),(label(3,j),j=1,number),dresid
      write(*,'(6d24.17)') dword,(label(1,j),j=1,number),
     | (label(2,j),j=1,number),(label(3,j),j=1,number),dresid

*     
*     write parameters values record (=a priori)
      read(20) dnum,(solveh(i),i=1,number)
      write(*,'(6d24.17)') dnum,(solveh(i),i=1,number)
    
*
*     write parameter a priori values record
      read(20) dnum,(solveh(i),i=1,number)
      write(*,'(6d24.17)') dnum,(solveh(i),i=1,number)

*
*     write normal equations
      do 90 i=1,number
      read(20) drow,delem,(norm(j*(j-1)/2+i),
     |  j=i,number),resid
      write(*,'(6d24.17)') drow,delem,(norm(j*(j-1)/2+i),
     |  j=i,number),resid
   90 continue
 
      end
