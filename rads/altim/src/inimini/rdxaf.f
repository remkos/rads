      subroutine rdxaf(iobs,rlat,rlon,it,cor,adj,ssh,func)
      implicit none
      include "data.inc"
      include "buffer.inc"
      integer ibuff,ios,readf,seekf,itrack,iobs,ip,icol,it,vpnt
      real    rlat,rlon,ssh_0,ssh,func(maxpar),cor,adj

      ip=(iobs-1)/bufsiz
      if (ip.ne.block) then
         block =ip
         minobs=block*bufsiz+1
         maxobs=min((block+1)*bufsiz,nrec)

         ios=seekf(fd,minobs*xaflen,0)
         ios=readf(fd,(maxobs-minobs+1)*xaflen,buffer)
      endif

      ibuff=(iobs-minobs)*xaflen+1
      call flxaf(buffer(ibuff),buffer(ibuff),rlat,rlon,
     |           itrack,ssh_0,func)

      it=itr(itrack)
      cor=0e0
      adj=0e0
      if (it.ne.0) then
         do icol=1,npar
            ip=vpnt(it,icol)
            adj=adj+adjst(ip)*func(icol)
            cor=cor+param(ip)*func(icol)
         enddo
      endif
      ssh=ssh_0-cor

      end

      subroutine flxaf(line2,line4,rlat,rlon,itrack,ssh_0,func)
      include "init.inc"
      include "data.inc"
      integer*2 line2(*)
      integer*4 line4(*),itrack,i
      real      rlat,rlon,ssh_0,func(maxpar),omega

      if (short) then
         rlat  =line4(2)/1e6
         rlon  =line4(3)/1e6
         ssh_0 =line4(4)/1e6
         itrack=wrap(int(line2(14)))
	 omega =(line4(6)-u0(itrack))*murad
         func(1)=1e0
	 func(2)=omega
	 do i=1,(npar-1)/2
            func(i*2  )=sin(i*omega)
            func(i*2+1)=cos(i*omega)
	 enddo
      else if (extend) then
         rlat  =line4(2)/1e6
         rlon  =line4(3)/1e6
         ssh_0 =line4(5)/1e6
         itrack=wrap(int(line2(14)))
         func(1)=1e0
         do i=2,npar
            func(i)=line4(6+i)/1e9
         enddo
      else
         rlat  =line4(2)/1e6
         rlon  =line4(3)/1e6
         ssh_0 =line4(5)/1e6
         itrack=wrap(int(line2(14)))
         do i=1,npar
            func(i)=line4(7+i)/1e9
         enddo
      endif

      end
