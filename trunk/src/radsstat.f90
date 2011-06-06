program radsstat

! This program reads the RADS data base and computes statistics
! by pass, cycle or day of a number of RADS data variables.
!
! usage: radsstat sat=<sat> [RADS_options] [options]
!-
! $Log: radsstat.f90,v $
! Created by Remko Scharroo, Altimetrics LLC, 2003/02/03
!-----------------------------------------------------------------------
use rads
use rads_time
use rads_misc
integer(fourbyteint) :: lstat=1
character(len=80) :: arg
character(len=32) :: wtype(0:3)=(/ &
	'box weight                  ', 'constant weight             ', &
	'area weighted               ', 'inclination-dependent weight'/)
character(len=512) :: format_string
integer(fourbyteint) :: nr=0, day=-999999, dayold=-99999, cycle, pass, i, j, k, l, per=0, wmode=0, nx, ny, kx, ky, ios
real(eightbytereal), allocatable :: z(:,:), lat_w(:)
real(eightbytereal) :: sini, step=1d0, x0, y0, dx=3d0, dy=1d0
type :: stat
	real(eightbytereal) :: wgt, mean, sum2, xmin, xmax
end type
type(stat), allocatable :: box(:,:,:), tot(:)
type(rads_sat) :: S
type(rads_pass) :: P

! Initialize RADS or issue help
if (iargc() < 1) call synopsis
call rads_init (S)
if (S%error /= rads_noerr) call synopsis

! If no sel= is given, use sla
if (S%nsel == 0)  call rads_parse_varlist (S, 'sla')

! Scan command line arguments
do i = 1,iargc()
	call getarg(i,arg)
	if (arg(:2) == '-d') then
		per = 0
		read (arg(3:),*,iostat=ios) step
	else if (arg(:2) == '-p') then
		per = 1
		read (arg(3:),*,iostat=ios) step
		step = dble(S%passes(2))/nint(S%passes(2)/step)
	else if (arg(:2) == '-c') then
		per = 2
		read (arg(3:),*,iostat=ios) step
	else if (arg(:2) == '-b') then
		wmode = 0
		call chartrans(arg(3:),'-+x',',,,')
		read (arg(5:),*,iostat=ios) dx,dy
	else if (arg(:2) == '-m') then
		wmode = 1
	else if (arg(:2) == '-a') then
		wmode = 2
	else if (arg(:2) == '-s') then
		wmode = 3
	else if (arg(:2) == '-l') then
		lstat = 2
	else if (arg(:4) == 'res=') then
		call chartrans(arg(5:),'-+x',',,,')
		read (arg(5:),*) dx,dy
	endif
enddo

! Set up the boxes
x0 = S%lon%info%limits(1)
y0 = S%lat%info%limits(1)
nx = int((S%lon%info%limits(2)-x0)/dx+0.999d0)
ny = int((S%lat%info%limits(2)-y0)/dy+0.999d0)
allocate (box(0:S%nsel,nx,ny),lat_w(ny))
allocate (tot(0:S%nsel))

! Set up the weights
sini = sin(S%inclination*rad)
if (wmode == 1) then
	lat_w = 1d0
else if (wmode == 3) then
	forall (ky=1:ny) lat_w(ky) = sqrt(max(1d-2,1d0-(sin((y0+(ky-0.5d0)*dy)*rad)/sini)**2))
else
	forall (ky=1:ny) lat_w(ky) = cos((y0+(ky-0.5d0)*dy)*rad)
endif

! Initialize statistics
box = stat(0d0, 0d0, 0d0, S%nan, S%nan)
call write_header

! Determine format for statistics
format_string = '(i9,f12.0'
do i = 1,S%nsel
	l = len_trim(format_string)
	! Add one decimal to the format
	arg = S%sel(i)%info%format(2:)
	call chartrans (arg,'.',' ')
	read (arg,*) j,k
	write (format_string(l+1:),'(",",i1,"(1x,f",i2.2,".",i2.2,")")') 2*lstat,j+1,k+1
enddo
l = len_trim(format_string)
format_string(l+1:) = ')'

! Start looping through cycles and passes
do cycle = S%cycles(1), S%cycles(2), S%cycles(3)
	! Process passes one-by-one
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cycle, pass)
		if (P%ndata > 0) call process_pass

		! Print the statistics at the end of the data pass (if requested)
		if (per == 1 .and. pass == nint(nint(pass/step)*step)) call print_stat
		call rads_close_pass (S, P)
	enddo

	! Print the statistics at the end of the cycle (if requested)
	if (per == 2 .and. cycle == nint(nint(cycle/step)*step)) call print_stat
enddo

! Flush the statistics (in "daily" mode) and close RADS4
call print_stat
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
call rads_this_is ('Revision: 4.0 $','Print RADS statistics per cycle, pass or day(s)')
call rads_synopsis ()
write (0,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -c[<n>]           : statistics per cycle or <n> cycles'/ &
'  -d[<n>]           : statistics per day (default) or <n> days'/ &
'  -p[<n>]           : statistics per pass or <n> passes'/ &
'  -b[dx,dy]         : average by boxes with size (default on: 3x1 degrees)'/ &
'  -m                : give all measurements equal weight'/ &
'  -a                : weight measurements by cosine of latitude'/ &
'  -s                : use inclination-dependent weight'/ &
'  -l                : print min and max in addition to mean and stddev'/ &
'  res=dx,dy         : size of averaging boxes (def: 3x1)')
stop
end subroutine synopsis

!***********************************************************************
! Process data for a single pass

subroutine process_pass
integer :: i, j

! Read the data for this pass
allocate (z(P%ndata,0:S%nsel))
z(:,0) = P%tll(:,1)
do j = 1,S%nsel
	call rads_get_var (S, P, S%sel(j), z(:,j))
enddo

! Update the statistics with data in this pass
do i = 1,P%ndata
	! Print the statistics (if in "daily" mode)
	day = floor(P%tll(i,1)/86400d0)
	if (per == 0 .and. day >= nint(dayold+step)) call print_stat
   	if (any(isnan(z(i,:)))) cycle ! Reject outliers

	! Update the box statistics
   	kx = floor((P%tll(i,2)-x0)/dx + 1d0)
   	ky = floor((P%tll(i,3)-y0)/dy + 1d0)
   	kx = max(1,min(kx,nx))
   	ky = max(1,min(ky,ny))
   	box(:,kx,ky)%wgt  = box(:,kx,ky)%wgt  + 1d0
   	box(:,kx,ky)%mean = box(:,kx,ky)%mean + z(i,:)
   	box(:,kx,ky)%sum2 = box(:,kx,ky)%sum2 + z(i,:)**2
   	box(:,kx,ky)%xmin = min(box(:,kx,ky)%xmin, z(i,:))
   	box(:,kx,ky)%xmax = max(box(:,kx,ky)%xmax, z(i,:))
   	nr = nr + 1
enddo

deallocate (z)
end subroutine process_pass

!***********************************************************************
! Print statistics for one batch of data

subroutine print_stat
integer(fourbyteint) :: j,yy,mm,dd
real(eightbytereal) :: w

if (nr == 0) then
	dayold = day
	return
endif

! Init global stats
tot = stat(0d0, 0d0, 0d0, S%nan, S%nan)

! Cycle trough all boxes and determine overall weighted mean
do ky=1,ny
	do kx=1,nx
		if (box(1,kx,ky)%wgt == 0d0) cycle ! Skip empty boxes

		! Determine the weight
		if (wmode == 0) then
	    	w = lat_w(ky) / box(1,kx,ky)%wgt
		else
	    	w = lat_w(ky)
		endif

		! Update overall statistics
		tot(:)%wgt  = tot(:)%wgt  + w * box(:,kx,ky)%wgt
		tot(:)%mean = tot(:)%mean + w * box(:,kx,ky)%mean
		tot(:)%sum2 = tot(:)%sum2 + w * box(:,kx,ky)%sum2
		tot(:)%xmin = min(tot(:)%xmin, box(:,kx,ky)%xmin)
		tot(:)%xmax = max(tot(:)%xmax, box(:,kx,ky)%xmax)
	enddo
enddo
	
! Divide by total weight to compute overall weighted mean and standard deviation

tot%mean = tot%mean / tot%wgt
if (nr > 1) then
	tot%sum2 = sqrt((tot%sum2 / tot%wgt - tot%mean**2) * nr / (nr - 1d0))
else
	tot%sum2 = S%nan
endif

! Print results
call mjd2ymd(dayold+46066,yy,mm,dd)
if (per == 0) then
	write (*,600) modulo(yy,100),mm,dd
else if (per == 1) then
	write (*,601) cycle,pass
else
	write (*,602) cycle,modulo(yy,100),mm,dd
endif
if (lstat == 1) then
	write (*,format_string) nr,tot(0)%mean,(tot(j)%mean,tot(j)%sum2,j=1,S%nsel)
else
	write (*,format_string) nr,tot(0)%mean,(tot(j)%mean,tot(j)%sum2,tot(j)%xmin,tot(j)%xmax,j=1,S%nsel)
endif

! Reset statistics
box = stat(0d0, 0d0, 0d0, S%nan, S%nan)
nr  = 0
dayold = day

600 format (3i2.2,$)
601 format (i3,i5,$)
602 format (i3,1x,3i2.2,$)
end subroutine print_stat

!***********************************************************************
! Write out the header

subroutine write_header
integer :: j0, j

600 format ('# Statistics of RADS variables (',a,')'/ &
'#'/'# Satellite : ',a,'/',a/'# Cycles    :',i5,' -',i5/ &
'# Passes    :',i5,' -',i5/'#'/'# Output columns:')
610 format ('#    ( 1) date [YYMMDD]')
611 format ('# ( 1, 2) cycle and pass')
612 format ('#    ( 1) cycle'/'#    ( 2) date at beginning of cycle [YYMMDD]')
620 format ('#    (',i2,') nr of measurements'/'#    (',i2,') mean time [',a,']')
621 format ('# (',i2,'-',i2,') mean and stddev of ',a,' [',a,']')
622 format ('# (',i2,'-',i2,') mean, stddev, min and max of ',a,' [',a,']')

write (*,600) trim(wtype(wmode)),trim(S%sat),trim(S%phase%name),S%cycles(1:2),S%passes(1:2)
if (per == 0) then
	write (*,610)
	j0 = 1
else if (per == 1) then
	write (*,611)
	j0 = 2
else
	write (*,612)
	j0 = 2
endif
write (*,620) j0+1,j0+2,trim(S%time%info%units)
do j = 1,S%nsel
	if (lstat == 1) then
		write (*,621) 2*j+j0+1,2*j+j0+2,trim(S%sel(j)%info%long_name),trim(S%sel(j)%info%units)
	else
		write (*,622) 4*j+j0+1,4*j+j0+2,trim(S%sel(j)%info%long_name),trim(S%sel(j)%info%units)
	endif
enddo
end subroutine write_header

!***********************************************************************

end program radsstat
