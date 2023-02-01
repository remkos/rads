program test_gottide

!  Test GOTTIDE

implicit double precision (a-h,o-z)
logical :: isdata

write(6,17)
time = 49046.803493d0
utc = (time-46066d0)*86400d0
dlat = -21.00
dlon = 5.0
call gotinit('GOT4.7',.false.)	! Init without load tide
do i=1,5
    call gottide (utc,dlat,dlon,tide,dummy)
    call lpetide (utc,dlat,1,tlp1,dummy)
    call lpeqmt (time*86400d0,dlat,tlp2)
    isdata=.not.isnan(tide)
    write(6,15) dlat,dlon,time,tide*1d2,tlp1*1d2,tlp2
    dlon = dlon + 5.d0
end do
dlon = 5.0
do i=1,8
    call gottide (utc,dlat,dlon,tide,dummy)
    call lpetide (utc,dlat,1,tlp1,dummy)
    call lpeqmt (time*86400d0,dlat,tlp2)
    isdata=.not.isnan(tide)
    write(6,15) dlat,dlon,time,tide*1d2,tlp1*1d2,tlp2
    time = time + 1.d0/24.d0
    utc = utc + 3600d0
end do
15 format(1x,2f10.2,f15.6,3f10.2)
17 format(5x,'latitude  longitude',7x,'time',6x,'tide    lptide'/)
end program test_gottide
