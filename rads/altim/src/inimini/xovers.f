      subroutine xovers
      implicit none
      include "data.inc"
      include "init.inc"
      include "stat.inc"
      include "buffer.inc"
      include "grid.inc"
      real    r,r2,rn,ssvar,ssvar2,edtmlt2,
     |        rlat,rlon,cor1,cor2,adj1,adj2,ssh1,ssh2,ts1,ts2
      real    func1(maxpar),func2(maxpar)
      integer ip,ios,openf,readf,closef,vpnt,mpnt,
     |        icol,irow,iobs,it1,it2,is1,is2,grdidx
      real    ej,e1,e2,wj,w1,w2
      character*4 ftype

* Open XXF file and initialise buffer maintenance

      if (fd.eq.0) then
	 if (inimode.le.0) then
	    fd=openf(xxf,'r+')
	 else
	    fd=openf(xxf,'r')
	 endif

         ios=readf(fd,4,ftype)
         ios=readf(fd,4,nrec)

         if (short) then
            if (ftype.ne.'@XXB') goto 1000
         else if (extend) then
            if (ftype.ne.'@XXE') goto 1000
         else 
            if (ftype.ne.'@XXF') goto 1000
         endif

	 bufsiz=nbytes/xxflen
	 block=-1
      endif

      edtmlt2=edtmlt**2
*     if (iter.eq.0) edtmlt2=1e30

* Read xovers one by one

      do iobs=1,nrec
         call rdxxf(iobs,rlat,rlon,it1,it2,cor1,cor2,
     |		adj1,adj2,ssh1,ssh2,func1,func2)

* Compute xover residual, sea surface variability, and edit criterion
* Reject areas with very high variability

         grdidx=nint((rlat-gy0)/gy)*gnx+nint((rlon-gx0)/gx)+1
         ssvar=grms(grdidx)
	 if (ssvar.gt.maxvar) ssvar=maxvar
	 if (it1.le.0 .or. it2.le.0 .or. inimode.le.0) goto 100
	 ssvar2=ssvar**2
	 ts1=sigma(it1)
	 ts2=sigma(it2)

         r=ssh1-ssh2
	 r2=r**2
	 rn=r2/edtmlt2

	 wj=2*ssvar2
	 ej=rn/wj

	 if (mpnt(it1,it2,1,1).ne.0 .and. ej.le.1.) then

* Joint residual, if both tracks are in the same matrix

            do icol=1,npar
               do irow=1,npar
                  ip=mpnt(it1,it2,icol,irow)
                  matrix(ip)=matrix(ip)-func1(icol)*func2(irow)/wj
               enddo
            enddo
	    e1=ej
	    e2=ej
	    w1=wj
	    w2=wj
	 else
	    w1=max(wj,ssvar2+ts2**2)
	    w2=max(wj,ssvar2+ts1**2)
	    e1=rn/w2
	    e2=rn/w1
	 endif

* Single residuals for matrix and vector (track 1)

	 if (e1.le.1.) then
            do icol=1,npar
               do irow=1,icol
                  ip=mpnt(it1,it1,irow,icol)
                  matrix(ip)=matrix(ip)+func1(icol)*func1(irow)/w1
               enddo
               ip=vpnt(it1,icol)
               vector(ip)=vector(ip)+r*func1(icol)/w1
            enddo
	 endif

* Single residuals for matrix and vector (track 2)

	 if (e2.le.1.) then
            do icol=1,npar
               do irow=1,icol
                  ip=mpnt(it2,it2,irow,icol)
                  matrix(ip)=matrix(ip)+func2(icol)*func2(irow)/w2
               enddo
               ip=vpnt(it2,icol)
               vector(ip)=vector(ip)-r*func2(icol)/w2
            enddo
	 endif

* Accumulate statistics for both single residuals

	    is1=satel(it1)
	    is2=satel(it2)
	    if (is1.le.is2) then
	       ip=h(is2)+is1
	    else
	       ip=h(is1)+is2
	    endif
            grdidx=nint((rlat-fy0)/fy)*fnx+nint((rlon-fx0)/fx)+1
	 if (e1.le.1.) then
	    s_rres(ip)=s_rres(ip)+r2
	    s_nres(ip)=s_nres(ip)+1
	    s_mcor(is1)=s_mcor(is1)+cor1
	    s_rcor(is1)=s_rcor(is1)+cor1**2
	    s_adjust(is1)=max(s_adjust(is1),abs(adj1))
	    s_ncor(is1)=s_ncor(is1)+1
            t_nres(it1)=t_nres(it1)+1
	    t_rres(it1)=t_rres(it1)+r2
	    trms(grdidx)=trms(grdidx)+r2
	    tnr(grdidx)=tnr(grdidx)+2
	 endif
	 if (e2.le.1.) then
	    s_rres(ip)=s_rres(ip)+r2
	    s_nres(ip)=s_nres(ip)+1
	    s_mcor(is2)=s_mcor(is2)+cor2
	    s_rcor(is2)=s_rcor(is2)+cor2**2
	    s_adjust(is2)=max(s_adjust(is2),abs(adj2))
	    s_ncor(is2)=s_ncor(is2)+1
            t_nres(it2)=t_nres(it2)+1
	    t_rres(it2)=t_rres(it2)+r2
	    trms(grdidx)=trms(grdidx)+r2
	    tnr(grdidx)=tnr(grdidx)+2
         endif

	 if (ej.lt.1.) then
	    t_eres(it1)=t_eres(it1)+1
	    t_eres(it2)=t_eres(it2)+1
	 endif
	 t_nrestot(it1)=t_nrestot(it1)+1
	 t_nrestot(it2)=t_nrestot(it2)+1

100	 continue
	 ts1=ssvar
	 ts2=ssvar
	 if (.not.good(it1)) ts1=-ts1
	 if (.not.good(it2)) ts2=-ts2
	 if (inimode.le.0) call wrxxf(iobs,ssh1,ssh2,ts1,ts2)

      enddo

      if (inimode.le.0) then
	 ios=closef(fd)
	 write (6,*) '... XXF file corrected'
	 fd=0
      endif
      return

 1000 call fin('inimini: file is not in XXB, XXE or XXF format')
      end
