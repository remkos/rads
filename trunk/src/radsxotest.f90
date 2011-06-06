!*radsxotest -- Test with RADS xover file
!+
program radsxotest
use typesizes
use netcdf
use rads_netcdf

character(len=80) :: filename
integer(fourbyteint) :: ncid, varid

call getarg(1,filename)

call nfs(nf90_open(filename, nf90_write, ncid))

call nfs(nf90_redef(ncid))
call nfs(nf90_def_var(ncid, 'satellite', nf90_int, (/ 1 /), varid))
call nfs(nf90_enddef(ncid))
call nfs(nf90_close(ncid))

end program radsxotest
