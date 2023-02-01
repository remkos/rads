* satcat.inc

* msat        maximum number of satellites
* name		satellite name
* ids		satellite international id
* mjd0,mjd1	mission period
* abias		altimeter bias
* orberr	orbit error estimate (sigma)
* maxerr	orbit error estimate (maximum)
* altsig	altimeter noise
* orbinc	orbital inclination
* orbalt	orbital altitude
* distmin	minimum distance between altimeter measurements

      integer	msat
      parameter (msat=20)
      character name(msat)*8
      integer	ids(msat),mjd0(msat),mjd1(msat)
      real  	abias(msat),orberr(msat),maxerr(msat),altsig(msat),
     |		orbinc(msat),orbalt(msat),distmin(msat)

      common /satcat/ name,ids,mjd0,mjd1,abias,orberr,maxerr,altsig,
     |	orbinc,orbalt,distmin
