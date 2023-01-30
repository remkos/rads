      subroutine gstide( dlat,dlon,time,tide,isdata )
*
*
*  name - gstide              name derivation - geosat derived tide
*
*  function -  to compute the ocean tidal height at a given time
*              and location.
*
*  language - fortran 77
*
*  arguments -
*     name      type  i/o               description
*     ----      ----  ---               -----------
*     dlat       d     i    north latitude (in degrees) for desired
*                           location.
*
*     dlon       d     i    east longitude (in degrees).
*
*     time       d     i    desired time, in seconds of mjd.  e.g.,
*                            jan 1, 1986 00:00 is 4011638400.
*
*     tide       d     o    computed tidal height, in cm.
*
*     isdata     l     o    logical denoting whether tide data exist at
*                           desired location.  if false, then tide is
*                           not modified.
*
*
*  usage notes -
*     using the input data as of nov 19, 1990, this routine computes
*     the height of the (ocean + load) tide at the desired location.
*     this quantity is appropriate for use in satellite altimetry.
*     the computed tide is composed of the 30 largest spectral lines
*     within both the diurnal and semidiurnal bands, sufficient for
*     representing all major constituents including nodal modulations.
*
*  processing logic -
*     the tidal height is computed at four grid points surrounding
*     the desired location, with the final result being a bilinear
*     interpolation.
*
*  file references -
*     two datasets are read on initial call to routine:
*     file 08 - is currently a set of 60 records denoting doodson
*              numbers for largest lines in diurnal & semidiurnal bands.
*     file 50 - reads global array of orthoweights which defines tidal
*              constants at all valid locations.
*
*  important local variables -
*     array uv is loaded with orthoweights as follows:
*         uv(i,j,l,m,n), where
*           i = 1 for u; 2 for v
*           j = 1,2 or 3  (3 complex coeffs per species.)
*           l = tidal species (1 or 2)
*           m = 1 to 244, longitude index number
*           n = 1 to 142, latitude index  (for lat -70.5 to +70.5).
*
*     programming note: this routine is written for general use, which
*     is not necessarily efficient for all applications.
*     if it is desired to compute tidal heights at, e.g., every 1-sec
*     observation along an arc, the program can be speeded up in
*     various ways... contact the authors for details.
*     also, no work has been done to vectorize this routine.
*     for machines like the cray, the double precision should be removed.
*
*  error processing - none.
*
*  technical references -
*     d. cartwright and r. ray, oceanic tides from geosat altimetry,
*        journal of geophysical research, 95, 3069-3090, 1990.
*        (particularly appendix a).
*
*  history -
*   version   date       programmer        change description
*   -------   ----       ----------        ------------------
*     1.0    8/17/90   ray/cartwright  initial version.
*     1.1   12/18/90   ray/cartwright  enhanced documentation;
*                                      use ascii input for portability.
*     1.2   12/31/91   r. ray        changed min w to 0.5 (old value of
*                                    0.2 allowed too much extrapolation)
*
*

      implicit double precision (a-h,o-z)
      real       uv(2,3,2,244,142), undef
      real       latmin,latmax,lonmin,lonmax
      dimension  u(3,2),v(3,2)
      save
      logical    init,isdata
      data       init/.true./, undef/99999./

*     on first call, read geosat-derived ortho-weights
*     ------------------------------------------------
      open (8,file='/u2/wi/tides/doodson.dat',form='formatted')
      open (50,file='/u2/wi/tides/orthowts.dat',form='formatted')
      if (init) then
         init = .false.
         dx = 360./244.
         do 4 j=1,142
         do 4 i=1,244
            uv(1,1,1,i,j) = undef
    4    continue
******   read(50) latmin,latmax,lonmin,lonmax        ! binary file
******   read(50) uv                                 ! not used.
         read(50,11) nx,ny,latmin,latmax,lonmin,lonmax
    8    continue
            read(50,12,end=10) i,j,slat,slon,u1,v1,u2,v2,u3,v3
            read(50,13)                      u4,v4,u5,v5,u6,v6
            uv(1,1,1,i,j) = u1
            uv(2,1,1,i,j) = v1
            uv(1,2,1,i,j) = u2
            uv(2,2,1,i,j) = v2
            uv(1,3,1,i,j) = u3
            uv(2,3,1,i,j) = v3
            uv(1,1,2,i,j) = u4
            uv(2,1,2,i,j) = v4
            uv(1,2,2,i,j) = u5
            uv(2,2,2,i,j) = v5
            uv(1,3,2,i,j) = u6
            uv(2,3,2,i,j) = v6
         go to 8
   10    continue
   11    format(2i4,4f8.3)
   12    format(2i4,2f8.3,6f8.3)
   13    format(24x,6f8.3)
      endif
      isdata = .true.

*     compute indices for desired position
*     ------------------------------------
      jlat1 = int(dlat - latmin) + 1
      jlat2 = jlat1 + 1
      if (jlat1.lt.1 .or. jlat2.gt.142) then
         isdata = .false.
         return
      endif
      xlon = dlon
      if (xlon.lt.lonmin) xlon = xlon + 360.d0
      ilon1 = int((xlon - lonmin)/dx) + 1
      ilon2 = ilon1 + 1
      if (ilon2.gt.244) ilon2 = 1
      if (ilon1.lt.1 .or. ilon1.gt.244) stop 301 ! should never happen.

      if (uv(1,1,1,ilon1,jlat1).eq.undef .and.
     *    uv(1,1,1,ilon2,jlat1).eq.undef .and.
     *    uv(1,1,1,ilon1,jlat2).eq.undef .and.
     *    uv(1,1,1,ilon2,jlat2).eq.undef)  then
         isdata = .false.
         return
      endif
      w1 = 0.d0
      w2 = 0.d0
      wx1 = (dx - (xlon - real(ilon1-1)*dx - lonmin))/dx
      wx2 = 1.d0 - wx1
      wy1 = 1.0 - dlat + real(jlat1-1) + latmin
      wy2 = 1.0 - wy1
*  interpolation weights:
*  w1,w2,w3,w4 are for northwest,northeast,southeast,southwest corners.
      w1 = wx1*wy2
      w2 = wx2*wy2
      w3 = wx2*wy1
      w4 = wx1*wy1
      s = 0.d0
      w = 0.d0

*     get orthoweights & compute tide
*     -------------------------------
      if (uv(1,1,1,ilon1,jlat1).ne.undef) then
         do 100 i=1,3
            u(i,1) = uv(1,i,1,ilon1,jlat1)
            v(i,1) = uv(2,i,1,ilon1,jlat1)
            u(i,2) = uv(1,i,2,ilon1,jlat1)
            v(i,2) = uv(2,i,2,ilon1,jlat1)
  100    continue
         call gstid1( u,v,time,tide1 )
         w = w4
         s = w4*tide1
      endif
      if (uv(1,1,1,ilon1,jlat2).ne.undef) then
         do 200 i=1,3
            u(i,1) = uv(1,i,1,ilon1,jlat2)
            v(i,1) = uv(2,i,1,ilon1,jlat2)
            u(i,2) = uv(1,i,2,ilon1,jlat2)
            v(i,2) = uv(2,i,2,ilon1,jlat2)
  200    continue
         call gstid1( u,v,time,tide1 )
         w = w + w1
         s = s + w1*tide1
      endif
      if (uv(1,1,1,ilon2,jlat2).ne.undef) then
         do 300 i=1,3
            u(i,1) = uv(1,i,1,ilon2,jlat2)
            v(i,1) = uv(2,i,1,ilon2,jlat2)
            u(i,2) = uv(1,i,2,ilon2,jlat2)
            v(i,2) = uv(2,i,2,ilon2,jlat2)
  300    continue
         call gstid1( u,v,time,tide1 )
         w = w + w2
         s = s + w2*tide1
      endif
      if (uv(1,1,1,ilon2,jlat1).ne.undef) then
         do 400 i=1,3
            u(i,1) = uv(1,i,1,ilon2,jlat1)
            v(i,1) = uv(2,i,1,ilon2,jlat1)
            u(i,2) = uv(1,i,2,ilon2,jlat1)
            v(i,2) = uv(2,i,2,ilon2,jlat1)
  400    continue
         call gstid1( u,v,time,tide1 )
         w = w + w3
         s = s + w3*tide1
      endif
      if (w.gt.0.5) then
         tide = s/w
      else
         isdata = .false.
      endif
      return
      end

      subroutine gstid1( u,v,ts,tide )

*  computes the tide at time ts from the orthoweights u,v.
*  this is a kernel routine, meant to be called only by gstide.

      implicit double precision (a-h,o-z)
      dimension   indx(30,2,4),amp(30,2),freq(30,2),
     *            pha(30,2),ct(30,2),st(30,2),at(3),bt(3),
     *            u00(2),u20(2),u21(2),u40(2),u41(2),v41(2),
     *            p(3,2),q(3,2),u(3,2),v(3,2),
     *            shpn(4),phc(4),dpd(4)
      logical     init
      save
*
      data phc/290.21d0, 280.12d0, 274.35d0, 343.51d0/,
     *     dpd/13.1763965d0,0.9856473d0,0.1114041d0,0.0529539d0/,
     *     u00/0.0298d0,0.0200d0/, u20/0.1408d0,0.0905d0/,
     *     u21/0.0805d0,0.0638d0/, u40/0.6002d0,0.3476d0/,
     *     u41/0.3025d0,0.1645d0/, v41/0.1517d0,0.0923d0/
      data tc/40431744.d2/, tslast/-9.99e10/
      data init/.true./
*
      if (init) then
      init = .false.
      rad = datan(1.d0)/45.d0
      dpld  = 36.d1 - dpd(1) + dpd(2)
      twopi = datan(1.d0)*8.d0
*
      fdpd  = 0.d0
      write(6,100)
  100 format(/'   semidiurnal & diurnal tidal potential terms:'/)
      do 201 m=1,2
         fdpd = fdpd + dpld
         do 101 l=1,30
            read(8,*) (indx(l,m,n),n=1,4), amp(l,m), pha(l,m)
            theta = fdpd
            do 1 n=1,4
               theta = theta + real(indx(l,m,n))*dpd(n)
    1       continue
            freq(l,m) = theta
            omegat = 2.d0*theta*rad
            ct(l,m) = 2.d0*dcos(omegat)
            st(l,m) = 2.d0*dsin(omegat)
            write(6,502) l,m, (indx(l,m,n),n=1,4), freq(l,m), amp(l,m),
     *                   pha(l,m)
  101    continue
         write(6,*) ' '
  201 continue
  502 format(1x,i3,i4,2x,4i2,f12.6,2(f10.5,f10.1,5x))
      endif

      if (ts.ne.tslast) then
         td = (ts - tc)/864.d2
*    compute 4 principal mean longitudes for a given time td
         do 107 n=1,4
            ph = phc(n) + td*dpd(n)
            shpn(n) = dmod(ph,36.d1)
  107    continue
         fd = td - dint(td)
         e  = 36.d1*fd - shpn(1) + shpn(2)
      endif
*
      otide = 0.d0
      do 300 m=1,2
         if (ts.ne.tslast) then
*
*  compute orthotides p(3),q(3) for given time ts in terms
*  of potential amplitudes at(3), bt(3)
            do 207 k=1,3
               at(k) = 0.d0
               bt(k) = 0.d0
  207       continue
            do 407 l=1,30
               ph = m*e + pha(l,m)
               do 307 n=1,4
                  i = indx(l,m,n)
                  if (i.eq.0) go to 307
                  ph = ph + real(i)*shpn(n)
  307          continue
               theta = ph*rad
               a = amp(l,m)*100.d0
               cr = a*dcos(theta)
               sr = a*dsin(theta)
               crct = cr*ct(l,m)
               crst = cr*st(l,m)
               srct = sr*ct(l,m)
               srst = sr*st(l,m)
               at(1) = at(1) + cr
               bt(1) = bt(1) - sr
               at(2) = at(2) + crct
               bt(2) = bt(2) - srct
               at(3) = at(3) + crst
               bt(3) = bt(3) - srst
  407       continue
            p(1,m) = u00(m)*at(1)
            q(1,m) = u00(m)*bt(1)
            p(2,m) = u20(m)*at(1) - u21(m)*at(2)
            q(2,m) = u20(m)*bt(1) - u21(m)*bt(2)
            p(3,m) = u40(m)*at(1) - u41(m)*at(2) - v41(m)*at(3)
            q(3,m) = u40(m)*bt(1) - u41(m)*bt(2) - v41(m)*bt(3)
         endif
         do 30 k=1,3
            otide = otide + u(k,m)*p(k,m) + v(k,m)*q(k,m)
   30    continue
  300 continue
      tide = otide
      tslast = ts
      return
      end
