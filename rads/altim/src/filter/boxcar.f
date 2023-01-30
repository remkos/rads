**BOXCAR -- Box car smoother
*+
      program boxcar

* This program performs a boxcar smoothing of a series of (t,z)
* data points. The two input columns are read from standard input.
* The smooth data are printed to standard output.
*
* syntax: boxcar [ length[s|y|d|n] ]
*
* The parameter 'length' determines the length of the boxcar.
* The length is given either is a number of records (add n),
* number of seconds (add s), number of days (add d), or a
* number of years (add y). Without length argument, no smoothing
* is applied.
*
* The input coordinate t is assumed to be is seconds. The
* additional s to the length parameter can be omitted.
*-
* $Log: boxcar.f,v $
* Revision 1.4  2008/02/18 02:27:06  rads
* - Changes for gfortran 4.3
*
* Revision 1.3  2005/03/31 17:38:47  remko
* - Make it work for decending series.
* - Skip lines starting with # or >
*
* Revision 1.2  2004/09/24 14:30:29  remko
* - Cosmetic changes only
*
* Revision 1.1  2004/09/21 21:05:16  remko
* - Filter routines added to altim/src suite
*
*-----------------------------------------------------------------------
      implicit none
      integer mrec,nrec/0/,i,j,n,l,lnblnk,nlength/999999/,iargc
      parameter (mrec=10000)
      real*8 a(mrec),b(mrec),suma,sumb,alength/1d30/
      character*80 arg

      do i=1,iargc()
         call getarg(i,arg)
         l=lnblnk(arg)
         if (arg(l:l).eq.'n') then
            read (arg(:l-1),*) nlength
         else if (arg(l:l).eq.'y') then
            read (arg(:l-1),*) alength
	    alength=alength*86400d0*365.25d0
         else if (arg(l:l).eq.'d') then
            read (arg(:l-1),*) alength
	    alength=alength*86400d0
         else if (arg(l:l).eq.'s') then
            read (arg(:l-1),*) alength
         else
            read (arg,*) alength
         endif
      enddo
      if (alength.eq.1d30 .and. nlength.eq.999999) nlength=1

10    read (*,550,end=100) arg
      if (arg(:1).eq.'#' .or. arg(:1).eq.'>') then
      	 ! skip
      else if (nrec.ge.mrec) then
	 call fin("Too many records")
      else
	 nrec=nrec+1
         read (arg,*) a(nrec),b(nrec)
      endif
      goto 10

100   nrec=nrec-1
      do i=1,nrec
         suma=0d0
	 sumb=0d0
	 n=0
	 do j=i,nrec
	    suma=suma+a(j)
	    sumb=sumb+b(j)
	    n=n+1
	    if (n.ge.nlength .or. abs(a(j)-a(i)).ge.alength) then
		write (*,*) suma/n,sumb/n
		goto 120
	    endif
	 enddo
120      continue
      enddo

550   format (a)
      end
