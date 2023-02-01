program testtide

! Driver for testing GOT and FES tide models and subroutines

!use typesizes
use tides
real(eightbytereal) :: lat,lon,utc,mjdd,otide,otide_sp,otide_lp=0d0,ltide,ltide_sp,ltide_lp=0d0,lpe,lpe_mf,ptide
integer(fourbyteint) :: mjd,i,j,mode=0,iargc
character(80) :: arg
logical :: air=.false.
type(gottideinfo) :: gotinfo
type(festideinfo) :: fesinfo
type(airtideinfo) :: airinfo
type(fes) :: fesinfo1, fesinfo2

mjd = 12053 + 33282
lon = -7.688d0
lat = 59.195d0

do i=1,iargc()
	call getarg(i,arg)
	if (arg(:4) == 'mjd=') then
		read (arg(5:),*) mjd
	else if (arg(:4) == 'lat=') then
		read (arg(5:),*) lat
	else if (arg(:4) == 'lon=') then
		read (arg(5:),*) lon
	else if (index(arg,'GOT') > 0) then
		mode = 1
		call gottideinit(arg,.true.,gotinfo)
	else if (index(arg,'FES2012/') > 0 .or. index(arg,'FES2014') > 0) then
		mode = 2
		j = fes_init(fesinfo1,fes_tide,fes_mem,arg)
		j = fes_init(fesinfo2,fes_radial,fes_mem,arg)
		call fes_set_nodal_time(fesinfo1,0d0)
		call fes_set_nodal_time(fesinfo2,0d0)
	else if (index(arg,'FES2012') > 0) then
		mode = 3
		call festideinit(arg,.false.,fesinfo)
	else if (index(arg,'FES') > 0) then
		mode = 3
		call festideinit(arg,.true.,fesinfo)
	else if (index(arg,'air') > 0) then
		call airtideinit(arg,airinfo)
		air = .true.
	endif
enddo

if (air) then
	write(*,10) 'JulDay','Hour','Latitude','Longitude','Short_tid','LP_tid','Pure_Tide','Geo_Tide','Rad_Tide','Pole_Tide','Air_Tide'
else
	write(*,10) 'JulDay','Hour','Latitude','Longitude','Short_tid','LP_tid','Pure_Tide','Geo_Tide','Rad_Tide','Pole_Tide'
endif
call poletideinit()
do i = 0,23
	utc = (mjd-46066)*86400d0+i*3600d0
	mjdd = mjd + i/24d0
	select case (mode)
	case (1)
		call gottide(gotinfo,utc,lat,lon,otide_sp,ltide_sp)
		call lpetide(utc,lat,0,otide_lp,lpe_mf)
	case (2)
		j = fes_eval(fesinfo1,utc,lat,lon,otide_sp,otide_lp)
		j = fes_eval(fesinfo2,utc,lat,lon,ltide_sp,ltide_lp)
	case (3)
		call festide(fesinfo,utc,lat,lon,otide_sp,otide_lp,ltide_sp,ltide_lp)
		call lpetide(utc,lat,1,lpe,lpe_mf)
		otide_lp = otide_lp + lpe - lpe_mf
	end select
	otide_sp = otide_sp * 100
	otide_lp = otide_lp * 100
	ltide_sp = ltide_sp * 100
	ltide_lp = ltide_lp * 100
	otide = otide_sp + otide_lp
	ltide = ltide_sp + ltide_lp
	ptide = poletide(mjdd,lat,lon)*1d2
	if (air) then
		write(*,20) (mjd-33282)+i/24d0,i,lat,lon,otide_sp,otide_lp,otide,otide+ltide,ltide,ptide,airtide(airinfo,utc,lat,lon)
	else
		write(*,20) (mjd-33282)+i/24d0,i,lat,lon,otide_sp,otide_lp,otide,otide+ltide,ltide,ptide
	endif
enddo

call gottidefree(gotinfo)
call festidefree(fesinfo)
call airtidefree(airinfo)
call poletidefree()

10 format(6x,a6,1x,a5,9(1x,a9))
20 format(f12.5,i6,9(1x,f9.3))
end program testtide
