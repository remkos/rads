* file: inclin.f
*
* INCLIN is the subroutine for inclination computation
*
      subroutine inclin(isat,isub)
      include "maxcom.inc"
      integer*4 isub,isub1,isub2,isat,i

      real*8 rlatc1,r1,xyz1(3)
      real*8 rlatc2,r2,xyz2(3)
      real*8 deltim,deltat,alpha,sina,u1,rnorm,scalar,arglat,deltau
      real*8 vector(3),incl,incl2

      isub1=1
      isub2=isub
      incl2=min(orbinc(isat),180-orbinc(isat))*rad/2

10    deltim=track(1,isub2)-track(1,isub1)
*
* Determine latitude in Geocentric system of first point: rlatc1
*
      call geocen(track(2,isub1),track(4,isub1)*1d3,rlatc1,r1)
      call polcar(rlatc1,track(3,isub1),r1,xyz1)
*
* Determine latitude in Geocentric system of second point: rlatc2,
* and also correct for the Earth rotation, t.i., converting to a
* non-rotating system
*
      call geocen(track(2,isub2),track(4,isub2)*1d3,rlatc2,r2)
      call polcar(rlatc2,track(3,isub2)+deltim*rotate,r2,xyz2)
*
* Compute scalar product of xyz1() and xyz2().
*
      scalar=xyz1(1)*xyz2(1)+xyz1(2)*xyz2(2)+xyz1(3)*xyz2(3)
      deltau=acos(scalar/r1/r2)
*
* When there is a significant time angle between the start and end point,
* try to reduce it.
* The reason for this is that if the angle between the first and the
* last point is in the order of 180 degrees it's very hard to compute
* proper node and inclination.
*
      i=isub2-isub1
      if (deltau.gt.0.6*pi .and. i.gt.2) then
	 if (abs(rlatc1).gt.incl2) isub1=isub1+1
	 if (abs(rlatc2).gt.incl2) isub2=isub2-1
	 if (isub2-isub1.ne.i) goto 10
      endif
*
* Compute vector normal on orbital plane
*
      vector(1)=xyz1(2)*xyz2(3)-xyz2(2)*xyz1(3)
      vector(2)=xyz1(3)*xyz2(1)-xyz2(3)*xyz1(1)
      vector(3)=xyz1(1)*xyz2(2)-xyz2(1)*xyz1(2)
*
* Compute angle between normal on equatorial plane (z=0 plane) and
* normal on orbital plane
*
      rnorm=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
      incl=acos(vector(3)/rnorm)
      progrd(itrknr)=(incl.lt.pi/2)
      if (.not.asc(itrknr)) incl=-incl
      inc(itrknr)=incl
*
* Compute
* - argument of latitude, measured from nearest node
* - angular velocity
* - time lap between beginning of track and node
* - longitude of the node
*
      theta(itrknr)=deltau/deltim
      u1=arglat(rlatc1,incl)
      deltat=-u1/theta(itrknr)
      sina=dtan(rlatc1)/dtan(incl)
      if (sina.gt.+1) then
         sina=+1
      else if (sina.lt.-1) then
         sina=-1
      endif
      alpha=dasin(sina)
      node(itrknr)=track(3,isub1)-alpha-deltat*rotate
      time(3,itrknr)=track(1,isub1)+deltat
      call geocen(track(2,1),track(4,1)*1d3,rlatc1,r1)
      call geocen(track(2,isub),track(4,isub)*1d3,rlatc2,r2)
      locat(1,itrknr)=track(3,1)
      locat(2,itrknr)=rlatc1
      locat(3,itrknr)=track(3,isub)
      locat(4,itrknr)=rlatc2
      u0(itrknr)=arglat(rlatc1,incl)
      end
