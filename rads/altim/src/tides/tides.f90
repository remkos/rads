module tides

use tides_air
use tides_fes
use tides_got
use tides_hret
use tides_lpe
use tides_ce
use tides_pole
use libfes

interface
	function webtideinit (dirname, ext, id)
	use typesizes
	character(*), intent(in) :: dirname, ext
	integer(fourbyteint) :: webtideinit, id
	end function webtideinit

	function webtide (id, utc, lat, lon, tide, tide2)
	use typesizes
	integer(fourbyteint) :: id
	real(eightbytereal), intent(in) :: utc, lat, lon
	real(eightbytereal), intent(out) :: tide, tide2
	integer(fourbyteint) :: webtide
	end function webtide

	subroutine webtidefree (id)
	use typesizes
	integer(fourbyteint) :: id
	end subroutine webtidefree
end interface

end module tides
