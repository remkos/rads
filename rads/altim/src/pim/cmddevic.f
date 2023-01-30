      subroutine cmddevic
      include "pim.inc"
      integer pgopen,lnblnk,l

* Open device

      if (popcmd('DEVIC',argum)) then
	 device=argum
	 call grtoup(argum,argum)
	 postscript=(index(argum,'/PS')+index(argum,'/VPS')+
     &          index(argum,'/CPS')+index(argum,'/VCPS').gt.0)
      endif
      l=lnblnk(device)
      write (0,'("Opening device ",a," ...")') device(:l)
      device_id=pgopen(device)
      if (device_id.le.0) stop "pim: can not open device"
      call pgsci(c_fg)
      call pgsch(def_ch)
      call pgslw(1)
      call pgsvp(0.,1.,0.,1.)

      end
