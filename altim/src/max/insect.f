**INSECT -- Guess location of xover between tracks J1 and J2
*
* INSECT calculates the estimates for the crossover locations using
* the ascending and descending nodes.
*-
* 16-Oct-1994 - Temp file is not written. Parse directly to XOVDET.
*
      subroutine insect(j1,j2)
      include "maxcom.inc"
      real*8 delinc,dellam
      integer*4 j1,j2,ja,jb,istep,maxstp
      real*8 rlatlo,rlathi,rlonlo,rlonhi,xdecl,xlam0,xlam1,
     |       ang1,ang2,dinc,decl,rlam,shiftlon
      parameter (delinc=1d-1,maxstp=5,dellam=1d-2*rad)
*
* Check inclination difference
*
      dinc=abs(inc(j1)-inc(j2))
      if (dinc.lt.delinc) goto 90
      dinc=abs(dinc-pi)
      if (dinc.lt.delinc) goto 90
*
* Determine whether area in common for both tracks exists
*
      rlatlo=max(min(locat(2,j1),locat(4,j1)),
     |           min(locat(2,j2),locat(4,j2)))
      rlathi=min(max(locat(2,j1),locat(4,j1)),
     |           max(locat(2,j2),locat(4,j2)))
      if (rlatlo.gt.rlathi) goto 90
      shiftlon=0
      rlonlo=max(min(locat(1,j1),locat(3,j1)),
     |           min(locat(1,j2),locat(3,j2)))
      rlonhi=min(max(locat(1,j1),locat(3,j1)),
     |           max(locat(1,j2),locat(3,j2)))
      if (rlonlo.le.rlonhi) goto 100
*
* If area is global, try again by shifting the node +360 or -360 degrees
*
      if (maxlon-minlon.lt.twopi) goto 90

      shiftlon=-twopi
      rlonlo=max(min(locat(1,j1),locat(3,j1))+shiftlon,
     |           min(locat(1,j2),locat(3,j2)))
      rlonhi=min(max(locat(1,j1),locat(3,j1))+shiftlon,
     |           max(locat(1,j2),locat(3,j2)))
      if (rlonlo.le.rlonhi) goto 100

      shiftlon=+twopi
      rlonlo=max(min(locat(1,j1),locat(3,j1))+shiftlon,
     |           min(locat(1,j2),locat(3,j2)))
      rlonhi=min(max(locat(1,j1),locat(3,j1))+shiftlon,
     |           max(locat(1,j2),locat(3,j2)))
      if (rlonlo.le.rlonhi) goto 100

90    continue
      nrxins=nrxins+1
      return
      
100   continue
*
* Shift the node by shiftlon
*
      node(j1)=node(j1)+shiftlon
*
* Determine track of satellite with largest inc for guessing algoritm.
* ang1 and ang2 are the angles between the tracks and the equator.
*
      ang1=abs(inc(j1))
      ang1=min(ang1,pi-ang1)
      ang2=abs(inc(j2))
      ang2=min(ang2,pi-ang2)
*
      if (ang1.gt.ang2) then
         ja=j1
         jb=j2
      else
         ja=j2
         jb=j1
      endif
*
* Check if inclinations of tracks are in same or opposit quadrants
* and start crossover guessing routine
*
      istep=0
      if ((progrd(j1).eqv.progrd(j2)).eqv.
     |       (asc(j1).eqv.asc(j2))) then
         xlam0=node(ja)
   20    xdecl=decl(jb,xlam0)
         xlam1=rlam(ja,xdecl)
         if (abs(xlam1-xlam0).gt.dellam .and. istep.lt.maxstp) then
            xlam0=xlam1
            istep=istep+1
            goto 20
         endif
*        write (6,1111) 'scharp: ',ja,jb,istep,
*    |                   inc(ja)/rad,inc(jb)/rad
*1111    format (a6,3i5,2(1x,f12.6))
      else
         xlam0=(node(ja)+node(jb))/2
   30    xdecl=decl(jb,xlam0)
         xlam1=(rlam(ja,xdecl)+xlam0)/2
         if (abs(xlam1-xlam0).gt.dellam .and. istep.lt.maxstp) then
            xlam0=xlam1
            istep=istep+1
            goto 30
         endif
*        write (6,1111) 'acute: ',ja,jb,istep,
*    |                   inc(ja)/rad,inc(jb)/rad
      endif
*
* XLAM1 and XDECL must be between the common area boundaries
*
*     write (6,'(2i4,2f10.3,2i4)') ja,jb,xlam0/rad,xlam1/rad,istep
*     if (xlam1.lt.rlonlo .or. xlam1.gt.rlonhi .or.
*    .    xdecl.lt.rlatlo .or. xdecl.gt.rlathi) return
*
* Shift the node back by shiftlon. Move crossover within area boundaries.
*
      node(j1)=node(j1)-shiftlon
      call grouplon(minlon,maxlon,xlam1)
*
* Now continue with XOVDET
*
      call xovdet(j1,j2,xdecl,xlam1)
      end
