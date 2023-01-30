!*MASKLIB -- A set of mask reading and query routines
!+
module masklib
use typesizes
type mask
    integer(fourbyteint) :: nx,ny
    real(eightbytereal) :: xmin,xmax,dx,ymin,ymax,dy
    integer (fourbyteint), allocatable :: mask(:)
end type
!
! Use the following routines to load masks into memory and query them:
! MASK_LOAD -- Load a mask into memory (BPM format)
! MASK_FREE -- Free mask buffer
! MASK_QUERY -- Look up value in buffered mask
! Information about the grid is stored in a structure of type 'mask'.
! Use the module 'masklib' to define the structure.
!-
! $Log: masklib.f90,v $
! Revision 1.1  2009/12/05 02:57:55  rads
! - Created by Remko Scharroo -- Altimetrics LLC
!
!-----------------------------------------------------------------------

contains

!&MASKBUFF -- Load a mask file into memory
!+
function mask_load (filenm, info)
integer :: mask_load
character(*), intent(in) :: filenm
type (mask), intent(inout) :: info
!
! This routine allocates memory and loads the contents of a mask
! file into the allocated memory.
!
! The mask file has to be coded in the PBM format. In the
! header of the PBM file the mask area can be coded using a line
! # AREA <lon0> <lon1> <lat0> <lat1>
! If this line is absent the default will be assumed:
! # AREA 0 360 -90 90
!
! All information about the mask is stored in the structure 'info'.
!
! After using mask_load, use mask_query to query the individual bits
! of the mask.
!
! Use mask_free to free the allocated memory. Note that this does not
! free the mask info structure.
!
! When the allocation of memory or the loading of the mask was
! unsuccessful, this will be reflected in the returned pointer value.
!
! Input argument:
!  filenm   : Name of the BPM file containing the mask.
!
! Output argument:
!  info     : Pointer to the mask structure
!
! Return value:
!  mask_load: 0 = No error
!             1 = Could not find or open mask
!             2 = Ilegal mask format
!             3 = Error loading mask
!-----------------------------------------------------------------------
character(256) :: line
integer(fourbyteint) :: fd,ios,nbytes,nwords,readf,openf
logical :: ltlend

! Open mask file. Check if file exists.

fd = openf(filenm,'r')
if (fd < 0) then
    write (0,1300) 'File not found:',trim(filenm)
    mask_load = 1
    return
endif

! Read header

ios = readf(fd,0,line) - 1
if (line(:ios) /= 'P4') then
    write (0,1300) 'Illegal mask format in',trim(filenm)
    mask_load = 1
    call closef(fd)
    return
endif

! Interpret mask header

info%xmin = 0d0
info%xmax = 360d0
info%ymin = -90d0
info%ymax = 90d0

do
    ios = readf(fd,0,line)-1
    if (ios <= 0) then
	exit
    else if (line(:6) == '# AREA') then
	read (line(7:ios),*) info%xmin,info%xmax,info%ymin,info%ymax
    else if (line(:1) /= '#') then
	read (line(:ios),*) info%nx,info%ny
	exit
    endif
enddo

! Determine dimensions and allocate memory

nbytes = (info%nx * info%ny + 7) / 8	! Total number of bytes of data in mask file
nwords = (nbytes + 3) / 4		! Determine number of 4-byte words (rounding up)
info%dx = (info%xmax - info%xmin) / info%nx
info%dy = (info%ymax - info%ymin) / info%ny
allocate (info%mask(0:nwords-1))

! Load entire data block.

ios = readf(fd,nbytes,info%mask)
call closef(fd)
if (ios /= nbytes) then
    write (0,1300) 'Error loading mask',trim(filenm)
    mask_load = 3
    return
endif

! Swap integers if machine is Little Endian.

if (ltlend()) call i4swap(nwords,info%mask(0:nwords-1))
mask_load = 0
return

! Error exits

1300 format ('mask_load: ',a,1x,a)

end function mask_load

!&MASK_FREE -- Free mask buffer
!+
subroutine mask_free (info)
type (mask), intent(inout) :: info
!
! This routine frees up the memory allocated by mask_load to store a mask.
! Note that this routine does not deallocate the mask structure 'info'
! itself, only the memory allocated to store the mask values.
!
! Input argument:
!  info : Mask info structure as returned by mask_load
!-----------------------------------------------------------------------
if (allocated(info%mask)) deallocate(info%mask)
end subroutine mask_free

!&MASK_QUERY -- Look up value in buffered mask
!+
function mask_query (info, x, y)
type(mask), intent(in) :: info
real(eightbytereal), intent(in) :: x,y
integer(fourbyteint) :: mask_query

! This function looks up the status of a single bit in a mask that was
! previously loaded using mask_load.
!
! The location at which the mask is interogated is given by x and y.
! These arguments are given in "world coordinates"; in other words, x
! must be between info%xmin and info%xmax, the coordinates of the left- and
! right-most mask point. Something similar holds for the y-coordinate.
!
! Upon exit, the function value mask_query will be either 1 or 0,
! depending on the value of the mask at the pixel in which (x, y)
! resides.
!
! When x and/or y are out of the limits of the mask a value -1 is
! returned.
!
! Input arguments:
!  info : Pointer to the mask structure as returned by GRIDBUFF
!  x, y : X- and Y-coordinate of the point to be interpolated
!
! Return value:
!  mask_query : Mask bit at the location (x, y)
!               0 = bit was not set
!               1 = bit was set
!              -1 = error occurred (e.g. x or y out of range)
!-----------------------------------------------------------------------
integer (fourbyteint) :: jx,jy,i

! If x or y are NaN or beyond the allowed range, return -1

if (x >= info%xmin .and. x <= info%xmax .and. y >= info%ymin .and. y <= info%ymax) then
else
    mask_query = -1
    return
endif
         
! Determine jx (0..nx-1) and jy (0..ny-1) of the closest bit cell

jx = min(info%nx - 1, int((x - info%xmin) / info%dx))
jy = min(info%ny - 1, int((info%ymax - y) / info%dy))

i = jy * info%nx + jx
if (btest(info%mask(i/32),31-mod(i,32))) then
    mask_query = 1
else
    mask_query = 0
endif
end function mask_query

end module masklib
