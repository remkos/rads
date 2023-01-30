**WRIXXF -- Write one xover record to file
*
* Write crossover to file in specified crossover XXF format if requested.
* Also ASCII file possible.
*-
* 16-Oct-1994 - Split off from XOVDET
*-----------------------------------------------------------------------
      subroutine wrixxf(rlatx,rlonx,ja,jb,utcxa,utcxb,
     .   seahxa,seahxb,orbhxa,orbhxb,auxa,auxb)
      real*8 rlatx,rlonx,utcxa,utcxb,seahxa,seahxb,orbhxa,orbhxb
      include "maxcom.inc"
      integer*4 ja,jb,i4xxf(12),i
      real*8 althxa,althxb,omegaa,omegab
      integer*2  it(2),isig(2),auxa(*),auxb(*)
      equivalence (i4xxf(5),it),(i4xxf(12),isig)

      if (specif(4:4).eq.'X') then
         omegaa=theta(ja)*(utcxa-time(3,ja))
         omegab=theta(jb)*(utcxb-time(3,jb))
         if (.not.asc(ja)) omegaa=omegaa+pi
         if (.not.asc(jb)) omegab=omegab+pi
         i4xxf(1)=nint(rlatx/murad)
         i4xxf(2)=nint(rlonx/murad)
         i4xxf(3)=nint(utcxa)
         i4xxf(4)=nint(utcxb)
         it(1)=ja
         it(2)=jb
         i4xxf(6)=nint(seahxa*1d6)
         i4xxf(7)=nint(seahxb*1d6)
         i4xxf(8)=i4xxf(6)
         i4xxf(9)=i4xxf(7)
         i4xxf(10)=nint(omegaa/murad)
         i4xxf(11)=nint(omegab/murad)
         isig(1)=nint(altsig(satid(ja))*1d3)
         isig(2)=nint(altsig(satid(jb))*1d3)
         write (wrunit,rec=nrxfnd+1) i4xxf
	 if (naux.gt.0)	write (wrunit+1,rec=nrxfnd+1)
     |		(auxa(i),auxb(i),i=1,naux)
      else if (specif(4:4).eq.'O') then
	 call wrxshort(wrunit,nrxfnd+1,rlatx/rad,rlonx/rad,
     |	    utcxa,utcxb,ja,jb,seahxa,seahxb,orbhxa,orbhxb,naux,auxa,auxb)
      else if (specif(4:4).eq.'S') then
	 call wrxshort(wrunit,nrxfnd+1,rlatx/rad,rlonx/rad,
     |	    utcxa,utcxb,ja,jb,seahxa,seahxb,-1d0,-1d0,0,0,0)
      else if (specif(4:4).eq.'A') then
	 althxa=orbhxa*1d3-seahxa
	 althxb=orbhxb*1d3-seahxb
	 call wrxasc(wrunit,utcxa,
     .			althxa,seahxa-seahxb,nrxfnd*2-1)
	 call wrxasc(wrunit,utcxb,
     .			althxb,seahxb-seahxa,nrxfnd*2)
      endif
      end

      subroutine wrxshort(iunit,irec,rlat,rlon,
     .utc1,utc2,j1,j2,h1,h2,orb1,orb2,naux,aux1,aux2)
      integer*4 iunit,irec,j1,j2,naux,i
      real*8 rlat,rlon,utc1,utc2,h1,h2,orb1,orb2
      integer*2 it(2),aux1(*),aux2(*)
      integer*4 i4xxs(9),i4xxo(2)
      equivalence (it,i4xxs(7))

      it(1)=j1
      it(2)=j2
      i4xxs(1)=nint(rlat*1d6)
      i4xxs(2)=nint(rlon*1d6)
      i4xxs(3)=int(utc1)
      i4xxs(4)=nint((utc1-i4xxs(3))*1d6)
      i4xxs(5)=int(utc2)
      i4xxs(6)=nint((utc2-i4xxs(5))*1d6)
      i4xxs(8)=nint(h1*1d6)
      i4xxs(9)=nint(h2*1d6)
      if (orb1.gt.0) then
         i4xxo(1)=nint(orb1*1d6)
         i4xxo(2)=nint(orb2*1d6)
         write (iunit,rec=irec) i4xxs,i4xxo
      else
         write (iunit,rec=irec) i4xxs
      endif
      if (naux.gt.0) write (iunit+1,rec=irec)
     |	 (aux1(i),aux2(i),i=1,naux)
      end

      subroutine wrxasc(iunit,utcx,althx,diffh,idod)
      integer*4 iunit,idod,mjd85,mjd
      real*8 utcx,althx,diffh,sec
      parameter (mjd85=46066)

      mjd=int(utcx/86400d0)
      sec=utcx-mjd*86400d0
      write (iunit,100) mjd+mjd85,sec,althx,diffh,idod
  100 format (i5,f13.6,f11.3,f8.3,i6)
      end
