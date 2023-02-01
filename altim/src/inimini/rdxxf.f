      subroutine rdxxf(iobs,rlat,rlon,it1,it2,cor1,cor2,
     |		adj1,adj2,ssh1,ssh2,func1,func2)
      implicit none
      include "data.inc"
      include "buffer.inc"
      integer ibuff,ios,readf,seekf,vpnt,
     |        iobs,itrack1,itrack2,it1,it2,ip,icol
      real    rlat,rlon,ssh1_0,ssh2_0,ssh1,ssh2,
     |        func1(maxpar),func2(maxpar),cor1,cor2,
     |	      adj1,adj2

      ip=(iobs-1)/bufsiz
      if (ip.ne.block) then
	 block =ip
	 minobs=block*bufsiz+1
         maxobs=min((block+1)*bufsiz,nrec)

         ios=seekf(fd,minobs*xxflen,0)
         ios=readf(fd,(maxobs-minobs+1)*xxflen,buffer)
      endif

      ibuff=(iobs-minobs)*xxflen+1
      call flxxf(buffer(ibuff),buffer(ibuff),rlat,rlon,
     |           itrack1,itrack2,ssh1_0,ssh2_0,func1,func2)

      it1=itr(itrack1)
      cor1=0e0
      adj1=0e0
      if (it1.ne.0) then
         do icol=1,npar
            ip=vpnt(it1,icol)
	    adj1=adj1+adjst(ip)*func1(icol)
            cor1=cor1+param(ip)*func1(icol)
         enddo
      endif
      ssh1=ssh1_0-cor1

      it2=itr(itrack2)
      cor2=0e0
      adj2=0e0
      if (it2.ne.0) then
         do icol=1,npar
            ip=vpnt(it2,icol)
	    adj2=adj2+adjst(ip)*func2(icol)
            cor2=cor2+param(ip)*func2(icol)
         enddo
      endif
      ssh2=ssh2_0-cor2

      end

      subroutine flxxf(line2,line4,rlat,rlon,
     |           itrack1,itrack2,ssh1_0,ssh2_0,func1,func2)
      include "init.inc"
      include "data.inc"
      integer*2 line2(*)
      integer line4(*)
      integer i,itrack1,itrack2
      real    rlat,rlon,ssh1_0,ssh2_0,omega1,omega2,
     |        func1(maxpar),func2(maxpar)

      if (short) then
         rlat   =line4(1)/1e6
         rlon   =line4(2)/1e6
         itrack1=wrap(int(line2(9)))
         itrack2=wrap(int(line2(10)))
         ssh1_0 =line4(6)/1e6
         ssh2_0 =line4(7)/1e6
         omega1 =(line4(10)-u0(itrack1))*murad
         omega2 =(line4(11)-u0(itrack2))*murad
	 func1(1)=1e0
	 func2(1)=1e0
	 func1(2)=omega1
	 func2(2)=omega2
	 do i=1,(npar-1)/2
            func1(i*2  )=sin(i*omega1)
            func1(i*2+1)=cos(i*omega1)
            func2(i*2  )=sin(i*omega2)
            func2(i*2+1)=cos(i*omega2)
         enddo
      else if (extend) then
         rlat =line4(1)/1e6
         rlon =line4(2)/1e6
         itrack1=wrap(int(line2(9)))
         itrack2=wrap(int(line2(10)))
         ssh1_0 =line4(6)/1e6
         ssh2_0 =line4(7)/1e6
         func1(1)=1e0
         func2(1)=1e0
         do i=2,npar
            func1(i)=line4(8+(i-1)*2)/1e9
            func2(i)=line4(9+(i-1)*2)/1e9
         enddo
      else
         rlat  =line4(1)/1e6
         rlon  =line4(2)/1e6
         itrack1=wrap(int(line2( 9)))
         itrack2=wrap(int(line2(10)))
         ssh1_0 =line4(6)/1e6
         ssh2_0 =line4(7)/1e6
         do i=1,npar
            func1(i)=line4(8+i*2)/1e9
            func2(i)=line4(9+i*2)/1e9
         enddo
      endif

      end
