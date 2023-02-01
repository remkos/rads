program hrettidetest
!
! This program tests the hrettide routine based on
! information provided by Ed Zaron
!
! Syntax: hrettidetest <pathname>
!
! where <pathname> is the pathname of the HRET NetCDF file
!-
! $Log: hrettidetest.f90,v $
! Revision 1.2  2019/03/08 09:57:29  rads
! - Introduced internal tide model HRET
!
! (c) Remko Scharroo - EUMETSAT
!-----------------------------------------------------------------------
use tides_hret
use typesizes

! Initialise

integer, parameter :: nt = 3
integer :: i
logical :: passed
character(len=512) :: pathname
type(hrettideinfo) :: info
type :: ref
	real(eightbytereal) :: time, lat, lon, tide_comp(6)
endtype
real(eightbytereal) :: time, tide, tide_comp(6)
type(ref) :: t(nt)

! Test cases

t(1) = ref(1046250908.335094d0, 18.948867d0, 149.491272d0, &
	(/-0.004648d0, -0.000328d0, 0.004537d0, 0.001753d0, 0.001196d0, 0.000073d0/))
t(2) = ref(1048578767.430392d0, 17.925634d0, 149.995549d0, &
	(/-0.001640d0, 0.001641d0, 0.000251d0, -0.000959d0, 0.000565d0, -0.001284d0/))
t(3) = ref(1051079248.994491d0, 18.960929d0, 149.837059d0, &
	(/-0.015131d0, -0.000893d0, -0.001623d0, 0.001529d0, 0.001067d0, -0.003253d0/))

! Load the model

call getarg(1, pathname)
call hrettideinit(pathname, info)

! Run the tests

do i = 1, nt
	time = (t(i)%time - 2446066.5d0) * 86400d0
	call hrettide(info, t(i)%time, t(i)%lat, t(i)%lon, tide, tide_comp)
	passed = all(abs(tide_comp - t(i)%tide_comp) < 1d-4)
	write (*,600) i, t(i)
	write (*,601) i, passed, tide_comp
enddo
600 format (i2,f20.6,2f12.6,6f10.6)
601 format (i2,l44,6f10.6)

! Deallocate memory

call hrettidefree(info)

end program hrettidetest
