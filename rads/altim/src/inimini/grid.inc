* grid.inc

* maxgrd	maximum grid size
* fnx,fny	output grid dimensions
* gnx,gny	input grid dimensions
* hnx,hny	input grid dimensions
* frms(maxgrd)	output variability grid
* grms(maxgrd)	input variability grid
* href(maxgrd)	reference surface grid
* trms(maxgrd)	temporary statistics grid
* tnr(maxgrd)	temporary statistics grid
* lon(2)	longitude boundaries
* lat(2)	latitude boundaries
* bound		boundaries (int*2)
* fx,fy		factors between grid points and degrees (output grid)
* gx,gy		factors between grid points and degrees (input grid)
* hx,hy		factors between grid points and degrees (ref grid)
* res		output grid name
* var		input grid name
* ref		reference surface grid name

      integer	maxgrd
      parameter (maxgrd=720*360)
      character*80	var,res,ref

      real  	grms(maxgrd),trms(maxgrd),href(maxgrd)
      integer*4 tnr(maxgrd),fnx,fny,gnx,gny,hnx,hny
      real	gx,gy,gx0,gx1,gy0,gy1,
     |		fx,fy,fx0,fx1,fy0,fy1,
     |		hx,hy,hx0,hx1,hy0,hy1
      integer*2 bound(4)

      real      lon(2),lat(2),minlon,maxlon,minlat,maxlat
      equivalence (minlon,lon(1)),(maxlon,lon(2))
      equivalence (minlat,lat(1)),(maxlat,lat(2))

      common /grid/  var,res,ref,
     |		     minlon,maxlon,minlat,maxlat,
     |		     gx,gy,gx0,gx1,gy0,gy1,gnx,gny,
     |		     fx,fy,fx0,fx1,fy0,fy1,fnx,fny,
     |		     hx,hy,hx0,hx1,hy0,hy1,hnx,hny,
     |		     grms,trms,href,tnr,bound
