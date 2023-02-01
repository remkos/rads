      subroutine surface
      implicit none
      include "stat.inc"
      include "init.inc"
      include "data.inc"
      include "grid.inc"
      include "buffer.inc"
      real    ts,ssvar,ssvar2,ssref,r,r2,cor,adj,
     |        rlat,rlon,ssh,func(maxpar)
      integer ios,openf,readf,closef,grdidx,vpnt,mpnt,
     |        it,is,ip,iobs,irow,icol
      real    e,w,gr1int4,edtmlt2
      character*4 ftype

      if (fd.eq.0) then
	 if (inimode.le.0) then
	    fd=openf(xaf,'r+')
	 else
	    fd=openf(xaf,'r')
	 endif

	 ios=readf(fd,4,ftype)
	 ios=readf(fd,4,nrec)

	 if (short) then
	    if (ftype.ne.'@XAB') goto 1000
	 else if (extend) then
	    if (ftype.ne.'@XAE') goto 1000
	 else 
	    if (ftype.ne.'@XAF') goto 1000
	 endif

	 bufsiz=nbytes/xaflen
	 block=-1
      endif

      edtmlt2=edtmlt**2
*     if (iter.eq.0) edtmlt2=1e30

      do iobs=1,nrec
         call rdxaf(iobs,rlat,rlon,it,cor,adj,ssh,func)

* Compute surface residual, sea surface variability, and edit criterion

         grdidx=nint((rlat-gy0)/gy)*gnx+nint((rlon-gx0)/gx)+1
         ssvar=min(maxvar,grms(grdidx))
	 if (ssvar.gt.maxvar) ssvar=maxvar
	 if (it.le.0 .or. inimode.le.0) goto 100
	 ssvar2=ssvar**2
	 ts=sigma(it)
	 w=ssvar2+ts**2

	 if (ref(1:1).eq.' ') then
	    r=ssh
	    r2=r**2
	    e=0
	 else
	    ssref=gr1int4(href,(rlon-hx0)/hx+1,(rlat-hy0)/hy+1,
     |		hnx,hny,hnx)
	    r=ssh-ssref
	    r2=r**2
	    e=r2/(edtmlt2*w)
	 endif

* Single residuals for matrix and vector

         if (e.le.1.) then
            do icol=1,npar
               do irow=1,icol
                  ip=mpnt(it,it,irow,icol)
                  matrix(ip)=matrix(ip)+func(icol)*func(irow)/w
               enddo
               ip=vpnt(it,icol)
               vector(ip)=vector(ip)+r*func(icol)/w
            enddo

* Accumulate statistics

            is=satel(it)
	    ip=h(is)+is
            s_rres(ip)=s_rres(ip)+r2
            s_nres(ip)=s_nres(ip)+1
            s_mcor(is)=s_mcor(is)+cor
            s_rcor(is)=s_rcor(is)+cor**2
            s_adjust(is)=max(s_adjust(is),abs(adj))
            s_ncor(is)=s_ncor(is)+1
            t_nres(it)=t_nres(it)+1
            t_rres(it)=t_rres(it)+r2
	    t_eres(it)=t_eres(it)+1
            grdidx=nint((rlat-fy0)/fy)*fnx+nint((rlon-fx0)/fx)+1
            trms(grdidx)=trms(grdidx)+r2
            tnr(grdidx)=tnr(grdidx)+2
	 endif

         t_nrestot(it)=t_nrestot(it)+1

* Write observation to file

100      continue
	 if (.not.good(it)) ssvar=-ssvar
	 if (inimode.le.0) call wrxaf(iobs,ssh,ssvar)

      enddo
      if (inimode.le.0) then
         ios=closef(fd)
	 write (6,*) '... XAF file corrected'
	 fd=0
      endif

      return

 1000 call fin('File is not in XAB, XAE or XAF format')
      end
