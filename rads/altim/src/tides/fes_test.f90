! This is a little program to test the Fortran interface to the FES library
!-----------------------------------------------------------------------

program fes_test
use typesizes
use libfes
use iso_c_binding

real(eightbytereal), parameter :: tide_ref(0:23) = (/ &
-100.990629d0, -137.104223d0, -138.482412d0, -104.345361d0, -42.515586d0, 32.374755d0, &
102.167627d0, 149.469495d0, 162.102726d0, 136.505573d0, 78.894879d0, 3.643044d0, &
-70.661218d0, -126.154394d0, -150.116619d0, -137.779885d0, -93.130691d0, -27.816052d0, &
41.546661d0, 97.256195d0, 124.945282d0, 117.469330d0, 77.027699d0, 14.659386d0 /)
real(eightbytereal), parameter :: lp_ref(0:23) = (/ &
0.903299d0, 0.876519d0, 0.849121d0, 0.821123d0, 0.792543d0, 0.763399d0, &
0.733711d0, 0.703496d0, 0.672774d0, 0.641566d0, 0.609890d0, 0.577768d0, &
0.545220d0, 0.512266d0, 0.478928d0, 0.445227d0, 0.411184d0, 0.376822d0, &
0.342161d0, 0.307223d0, 0.272032d0, 0.236608d0, 0.200974d0, 0.165154d0 /)
real(eightbytereal), parameter :: load_ref(0:23) = (/ &
3.881161d0, 4.328335d0, 3.710694d0, 2.134257d0, -0.052047d0, -2.341404d0, &
-4.194242d0, -5.171971d0, -5.045669d0, -3.852375d0, -1.884925d0, 0.381964d0, &
2.410565d0, 3.733913d0, 4.070741d0, 3.392764d0, 1.927624d0, 0.097163d0, &
-1.592943d0, -2.683946d0, -2.881870d0, -2.132584d0, -0.635177d0, 1.209534d0 /)
real(eightbytereal), parameter :: lon = -7.688d0, lat = 59.195d0
character(len=80) :: ini
integer(fourbyteint) :: total_err = 0

call getarg (1, ini)
call test(fes_io)
call test(fes_mem)
write (*,*) 'Total number of errors: ', total_err
if (total_err /= 0) stop

contains

subroutine test (access)
integer(fourbyteint), intent(in) :: access
integer(fourbyteint) :: err, hour
real(eightbytereal) :: tide, time, lp, load, loadlp
type(fes) :: short_tide, radial_tide

600 format ('*** testing libfes with ',a,'...')
if (access == fes_io) then
	write (*,600) "direct access"
else if (access == fes_mem) then
	write (*,600) "memory access"
else
	write (*,600) "buffered access"
endif

write (*,*) 'Initialise ocean tides'
err = fes_new (short_tide, fes_tide, access, trim(ini) // c_null_char)
call fes_perror (short_tide)
write (*,*) 'Initialise load tides'
err = fes_new (radial_tide, fes_radial, access, trim(ini) // c_null_char)
call fes_perror (radial_tide)

do hour = 0, 23
	time = 12053d0 + hour / 24d0
	err = fes_core (short_tide, lat, lon, time, tide, lp)
	call fes_perror (short_tide)
	call check_tide (tide, tide_ref, hour, 'tide')
	call check_tide (lp, lp_ref, hour, 'lp')
	err = fes_core (radial_tide, lat, lon, time, load, loadlp)
	call fes_perror (radial_tide)
	call check_tide (load, load_ref, hour, 'load')
enddo
end subroutine test

subroutine check_tide (val, ref, hour, txt)
real(eightbytereal), intent(in) :: val, ref(0:)
integer(fourbyteint), intent(in) :: hour
character(len=*), intent(in) :: txt
if (abs(val - ref(hour)) > 1d-5) then
	write (*,*) 'hour = ', hour, ', type = ', txt, ': ', val, ref(hour)
	total_err = total_err + 1
endif
end subroutine check_tide

end program fes_test
