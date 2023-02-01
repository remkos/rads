      subroutine iniinit
      implicit none
      include "init.inc"
      include "grid.inc"
      include "data.inc"
      include "satcat.inc"
      integer	i,l,gridrd4
      real	z0,z1

      character*128 arg,error/' '/
      integer	iargc,iarg,l0,l1,l2,l3,l4,l5,l6

      namelist /satcat_nml/ name,ids,mjd0,mjd1,abias,orberr,maxerr,
     |	altsig,orbinc,orbalt,distmin
      namelist /inimini_nml/ edtmlt,perc,maxiter,bndwth,edttrk,
     |	xtf,xxf,xaf,lon,lat,fx,var,res,ref,init,smcells,maxvar

* Initialises some variables at the start of the run

* Read satcat.nml

      arg='/user/altim'
      call checkenv('ALTIM',arg,l)
      arg(l+1:)='/nml/satcat.nml'
      open (9,file=arg,status='old',err=10)
      read (9,satcat_nml)
      close (9)
10    continue
      open (9,file='satcat.nml',status='old',err=11)
      read (9,satcat_nml)
      close (9)
11    continue

* Read inimini.nml

      arg(l+1:)='/nml/inimini.nml'
      open (9,file=arg,status='old',err=20)
      read (9,inimini_nml)
      close (9)
20    continue
      open (9,file='inimini.nml',status='old',err=21)
      read (9,inimini_nml)
      close (9)
21    continue

* If any of the file names start with $ALTIM, replace with the environment variable

      call altim_strip(arg(:l),l,xaf)
      call altim_strip(arg(:l),l,xxf)
      call altim_strip(arg(:l),l,xtf)
      call altim_strip(arg(:l),l,res)
      call altim_strip(arg(:l),l,ref)
      call altim_strip(arg(:l),l,var)

* Read command line arguments

      do iarg=1,iargc()
	 call getarg(iarg,arg)
	 if (arg(1:7).eq.'edtmlt=') then
	    read (arg(8:),*) edtmlt
	 else if (arg(1:5).eq.'perc=') then
	    read (arg(6:),*) perc
	 else if (arg(1:8).eq.'maxiter=') then
	    read (arg(9:),*) maxiter
	 else if (arg(1:7).eq.'bndwth=') then
	    read (arg(8:),*) bndwth
	 else if (arg(1:4).eq.'lon=') then
	    read (arg(5:),*) lon
	 else if (arg(1:4).eq.'lat=') then
	    read (arg(5:),*) lat
	 else if (arg(1:3).eq.'fx=') then
	    read (arg(4:),*) fx
	 else if (arg(1:4).eq.'xaf=') then
	    xaf=arg(5:)
	 else if (arg(1:4).eq.'xxf=') then
	    xxf=arg(5:)
	 else if (arg(1:4).eq.'xtf=') then
	    xtf=arg(5:)
	 else if (arg(1:4).eq.'res=') then
	    res=arg(5:)
	 else if (arg(1:4).eq.'var=') then
	    var=arg(5:)
	 else if (arg(1:4).eq.'ref=') then
	    ref=arg(5:)
	 else if (arg(1:7).eq.'maxvar=') then
	    read (arg(8:),*) maxvar
	 else if (arg(1:8).eq.'smcells=') then
	    read (arg(9:),*) smcells
	 else if (arg(1:5).eq.'init=') then
	    read (arg(6:),*) init
	 else if (arg(1:7).eq.'edttrk=') then
	    read (arg(8:),*) edttrk
	 else
	    l0=index(arg,' ')-1
	    xaf=arg(:l0)//'.xaf'
	    xtf=arg(:l0)//'.xtf'
	    xxf=arg(:l0)//'.xxf'
	 endif
      enddo

* Check calling program

      call getarg(0,arg)
      if (index(arg,'iniedit').gt.0) then
         inimode=0
      else if (index(arg,'inimini').gt.0) then
	 inimode=1
	 xaf=' '
      else if (index(arg,'inisurf').gt.0) then
	 inimode=2
	 xxf=' '
      else
         call fin('unknown call to inimini')
      endif

* Check definition of files

      if (xtf.eq.' ' .or. (xxf.eq.' ' .and. xaf.eq.' ')) then
         error='inimini: specify xtf and xaf/xxf'
	 goto 900
      endif

* Get sizes

      call get_size

      if (maxlon.eq.minlon .or. maxlat.eq.minlat) then
	 error='inimini: specify lat and lon'
	 goto 900
      endif

* Some math

      pi=4*atan(1e0)
      rad=pi/180e0
      murad=rad/1e6
      sqrt_two=sqrt(2e0)

* Set wrap function

      do i=1,32767
         wrap(i)=i
      enddo
      do i=-32768,0
         wrap(i)=65536+i
      enddo

* Get the global xover rms grid

      gnx=0
      gny=maxgrd
      if (var.eq.'-') then
        gnx=2
        gny=2
        gx0=-180
        gx1=180
        gy0=-90
        gy1=90
      else
        i=gridrd4(var,gnx,gny,grms,gx0,gx1,gy0,gy1,z0,z1)
        if(i.ne.0)call fin('inimini: error loading variability grid')
      endif
      gx=(gx1-gx0)/(gnx-1)
      gy=(gy1-gy0)/(gny-1)
      do i=1,gnx*gny
	 grms(i)=max(0.01,grms(i))
      enddo

* Get sea surface height reference grid (only when specified)

      if (inimode.ne.2) ref=' '
      if (ref(1:1).ne.' ') then
         hnx=0
         hny=maxgrd
         i=gridrd4(ref,hnx,hny,href,hx0,hx1,hy0,hy1,z0,z1)
         if (i.ne.0)call fin('inimini: error loading surface grid')
         hx=(hx1-hx0)/(hnx-1)
         hy=(hy1-hy0)/(hny-1)
      endif

* Set the sizes for the output grid

      fy=fx
      fnx=nint((maxlon-minlon)/fx)
      fny=nint((maxlat-minlat)/fy)
      fx=(maxlon-minlon)/fnx
      fy=(maxlat-minlat)/fny
      fx0=minlon+fx/2
      fx1=maxlon-fx/2
      fy0=minlat+fy/2
      fy1=maxlat-fy/2

*     call errset(2104,20,20,0,0,0)

* Write inputs

900   l0=index(arg,' ')-1
      l1=index(xtf,' ')-1
      l2=index(xxf,' ')-1
      l3=index(xaf,' ')-1
      l4=index(res,' ')-1
      l5=index(var,' ')-1
      l6=index(ref,' ')-1

      write (6,1010) arg(:l0),version,
     |	minlon,maxlon,minlat,maxlat,var(:l5),maxvar,init,xtf(:l1)
      if (inimode.ne.2) write (6,1020) xxf(:l2)
      if (inimode.ne.1) write (6,1030) xaf(:l3)
      if (inimode.eq.2) write (6,1035) ref(:l6)
      if (inimode.ne.0) write (6,1040)
     |	edtmlt,edttrk,bndwth,res(:l4),fx,smcells,
     |  perc,maxiter

1010  format (
     |'This is ',a,', version ',f6.1//
     |'Options (command line or inimini.nml or',
     |' $ALTIM/nml/inimini.nml) : '/
     |'Longitude                 :      lon=',2f9.3/
     |'Latitude                  :      lat=',2f9.3/
     |'Input variability grid    :      var=',a/
     |'Maximum variability       :   maxvar=',f9.3/
     |'Initialise parameters     :     init=',l9/
     |'XTF filename              :      xtf=',a)
1020  format (
     |'XXF filename              :      xxf=',a)
1030  format (
     |'XAF filename              :      xaf=',a)
1035  format (
     |'Reference surface grid    :      ref=',a)
1040  format (
     |'Residual edit multiplier  :   edtmlt=',f9.3/
     |'Track edit multiplier     :   edttrk=',f9.3/
     |'Max bandmatrix halfwidth  :   bndwth=',i9/
     |'Output variability grid   :      res=',a/
     |'Degrees per RMS cell      :       fx=',f9.3/
     |'Nr cells for smoothing    :  smcells=',i9/
     |'Percentage criterium      :     perc=',f9.3/
     |'Maximum iterations        :  maxiter=',i9)

      if (error.ne.' ') call fin(error)
      end

      subroutine altim_strip (altim,l,arg)
      character*(*) altim,arg
      integer l
      if (arg(:6).eq.'$ALTIM') arg=altim(:l)//arg(7:)
      end
