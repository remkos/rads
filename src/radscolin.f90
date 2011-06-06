!*radscolin -- Make collinear data sets from RADS
!+
program radscolin
!
! This program provides a quick and crude way to make collinear data
! sets (in ASCII).
!
! Usage: radscolin [RADS_options] [options]
!-
! $Log: radscolin.f90,v $
! Created by Remko Scharroo, Altimetrics LLC, 2002-01-11
!-----------------------------------------------------------------------
use rads
integer(fourbyteint), parameter :: msat = 10
type(rads_sat) :: S(msat)
integer(fourbyteint) :: nsel = 0, reject = 9999, cycle, pass, i, j, ios, &
	nbins, nsat = 0, ntrx = 0, ntrx0, ntrx1, ptrx0, ptrx1, type_sla = 1, step = 1
real(eightbytereal) :: dt = 0.97d0
character(len=80) :: arg
character(len=512) :: format_string
logical :: numbered = .false., counter = .false.
real(eightbytereal), allocatable :: data(:,:,:)
integer(fourbyteint), allocatable :: nr_in_bin(:)
type :: stat
	integer(fourbyteint) :: nr
	real(eightbytereal) :: mean,sum2
end type
type(stat), allocatable :: pp(:,:)

! Initialize RADS or issue help
if (iargc() < 1) call synopsis
S%sat = '' ! Initialize blank
S%error = rads_noerr
call rads_init (S)
if (any(S%error /= rads_noerr)) call synopsis

! If no sat= is given, exit
if (S(1)%sat == '') call rads_exit ('Need at least one sat= option')

! Determine how many satellites and cycles.
! Also check that the same number of variables have been selected for each satellite
do i = 1,msat
	if (S(i)%sat == '') exit
	if (S(i)%nsel == 0) call rads_parse_varlist (S(i), 'sla')
	if (S(i)%nsel /= S(1)%nsel) call rads_exit ('Unequal amount of variables on sel= for different satellites')
	if (any(S(i)%passes /= S(1)%passes)) call rads_exit ('Unequal number of passes per cycle for different satellites')
	ntrx = ntrx + (S(i)%cycles(2) - S(i)%cycles(1)) / S(i)%cycles(3) + 1
	dt = max(dt,S(i)%dt1hz)
	nsat = i
enddo
nsel = S(1)%nsel

! Determine which column to check for NaN (default: 1st)
do j = 1,nsel
	if (S(nsat)%sel(j)%info%datatype == rads_type_sla) type_sla = j
enddo

! Set default column ranges
ntrx0 = 1
ntrx1 = ntrx
ptrx0 = 1
ptrx1 = ntrx
reject = ntrx

! Scan command line arguments
do i = 1,iargc()
	call getarg(i,arg)
	if (arg(:3) == '-rn') then
		reject = ntrx
	else if (arg(:2) == '-r') then
		reject = 0
		read (arg(3:),*,iostat=ios) reject
	else if (arg(:5) == 'step=') then
		read (arg(6:),*) step
	else if (arg(:3) == 'dt=') then
		read (arg(4:),*) dt
	else if (arg(:2) == '-a') then
		ptrx1 = ntrx + 1
	else if (arg(:2) == '-A') then
		ptrx0 = ntrx + 1
		ptrx1 = ntrx + 1
	else if (arg(:2) == '-s') then
		ptrx1 = ntrx + 2
	else if (arg(:2) == '-S') then
		ptrx0 = ntrx + 1
		ptrx1 = ntrx + 2
	else if (arg(:2) == '-n') then
		numbered = .true.
	else if (arg(:2) == '-N') then
		counter = .true.
	endif
enddo
reject = max(0,min(reject,ntrx))

! Build format string
if (ptrx1 > ptrx0) then
	write (format_string,'("(",a,",",i3,"(1x,",a,"),")') &
		trim(S(nsat)%sel(1)%info%format),(ptrx1-ptrx0),trim(S(nsat)%sel(1)%info%format)
else
	write (format_string,'("(",a,",")') trim(S(nsat)%sel(1)%info%format)
endif
do i = 2,nsel
	write (format_string(len_trim(format_string)+1:),'(i3,"(1x,",a,"),")') &
		(ptrx1-ptrx0+1),trim(S(nsat)%sel(i)%info%format)
enddo
format_string(len_trim(format_string):) = ')'

! Allocate data arrays
nbins = nint(S(1)%phase%repeat_days/S(1)%phase%repeat_passes*86400/dt + 60d0) ! Number of bins
allocate (data(ntrx+2,nsel,-nbins/2:nbins/2), nr_in_bin(-nbins/2:nbins/2), pp(ptrx0:ptrx1,nsel))

! Read one pass for each satellites at a time
do pass = S(1)%passes(1), S(1)%passes(2), S(1)%passes(3)
	call process_pass
enddo

contains

!***********************************************************************

subroutine synopsis
call rads_this_is ('$Revision: 4.0 $','Make collinear data sets from RADS')
call rads_synopsis ()
write (0,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r#               : reject data when there are fewer than # tracks with valid SLA'/ &
'                      (default: # = number of selected cycles)'/ &
'  -r0, -r           : keep all stacked data points'/ &
'  -rn               : reject data when any track is NaN (default)'/ &
'  dt=dt             : set minimum bin size in seconds (default is determined by satellite)'/ &
'  -a                : print mean in addition to pass data'/ &
'  -A                : print only mean'/ &
'  -s                : print mean and standard deviation in addition to pass data'/ &
'  -S                : print only mean and standard deviation'/ &
'  -n                : add record number' / &
'  -N                : add number of measurements in each bin')
end subroutine synopsis

!***********************************************************************
! Process the data for a single pass

subroutine process_pass
real(eightbytereal), allocatable :: temp(:)
integer, allocatable :: bin(:)
integer :: i, j, k, l
type(rads_pass) :: P
type(stat) :: pm

! Initialize data array
data = S(1)%nan
i = 0
do j = 1,nsat
	do cycle = S(j)%cycles(1), S(j)%cycles(2), S(j)%cycles(3)
		i = i + 1 ! track id
		call rads_open_pass (S(j), P, cycle, pass)
		if (P%ndata > 0) then
			allocate (temp(P%ndata),bin(P%ndata))
			bin = nint((P%tll(:,1) - P%equator_time) / dt) ! Store bin nr associated with measurement
			do l = 1,nsel
				call rads_get_var (S(j), P, S(j)%sel(l), temp)
				data(i,l,bin(:)) = temp(:)
			enddo
			deallocate (temp,bin)
		endif
		call rads_close_pass (S(j), P)
	enddo
enddo

! Check if there are common bins
nr_in_bin = 0
do k = -nbins/2,nbins/2,step
	j = count(.not.isnan(data(:,type_sla,k)))
	if (j >= reject) nr_in_bin(k) = j
enddo
if (sum(nr_in_bin) == 0) return

! Print the header
call write_header

! Do per-measurement and per-pass stats
call begin_stat (pp)
do k = -nbins/2,nbins/2,step
	if (nr_in_bin(k) == 0) cycle
	do l = 1,nsel
		call begin_stat (pm)
		do i = ntrx0,ntrx1
			call update_stat (pm, data(i,l,k))
		enddo
		call end_stat (pm)
		data (ntrx+1,l,k) = pm%mean
		data (ntrx+2,l,k) = pm%sum2
		call update_stat (pp(ptrx0:ptrx1,l), data(ptrx0:ptrx1,l,k))
	enddo
enddo
call end_stat (pp)

! Print out data that are common to some passes
635 format(i6)

do k = -nbins/2,nbins/2,step
	if (nr_in_bin(k) == 0) cycle
	write (*,format_string,advance='no') data(ptrx0:ptrx1,:,k)
	if (counter) write (*,635,advance='no') nr_in_bin(k)
	if (numbered) write (*,635,advance='no') k
	write (*,*)
enddo

! Write per-pass stats
640 format('# ',a,': ')
645 format(i4,i5)
650 format('# nr : ',400i6)

write (*,640,advance='no') 'avg'
write (*,format_string,advance='no') pp%mean
write (*,645) S(1)%cycles(1),pass
write (*,640,advance='no') 'std'
write (*,format_string,advance='no') pp%sum2
write (*,645) S(1)%cycles(1),pass
write (*,650) pp%nr,S(1)%cycles(1),pass

end subroutine process_pass

!***********************************************************************
! Generate statistics

elemental subroutine begin_stat (s)
type(stat), intent(inout) :: s
s = stat (0, 0d0, 0d0)
end subroutine begin_stat

elemental subroutine update_stat (s, x)
type(stat), intent(inout) :: s
real(eightbytereal), intent(in) :: x
real(eightbytereal) :: q, r
if (isnan(x)) return
s%nr = s%nr + 1
q = x - s%mean
r = q / s%nr
s%mean = s%mean + r
s%sum2 = s%sum2 + r * q * (s%nr - 1)
end subroutine update_stat

elemental subroutine end_stat (s)
type(stat), intent(inout) :: s
if (s%nr == 0) s%mean = s%mean / s%mean ! To make NaN
s%sum2 = sqrt(s%sum2/(s%nr-1))
if (s%nr <= 1) s%sum2 = s%sum2 / s%sum2 ! To make NaN
end subroutine end_stat

!***********************************************************************
! Write the pass header

subroutine write_header
logical :: first = .true.

605 format('# Satellite data selections:')
610 format('#   sat=',a,1x,'cycle=',i3.3,'-',i3.3,' pass=',i4.4)
611 format(' (',i3,'-',i3,')')
615 format('#   ',a,' (',i3,')')
620 format('#'/'# Column ranges for each variable:')
622 format('# ',i4,' -',i4,' : ',a,' [',a,']')
625 format('# ',i4,7x,': ',a)
630 format('#')

if (.not.first) write (*,*) ! Skip line between passes
first = .false.

! Describe data set per variable
write (*,605)
i = 1
do j = 1,nsat
	if (ptrx0 == 1) then !
		write (*,610,advance='no') S(j)%sat,S(j)%cycles(1:2),pass
		write (*,611) i,i + (S(j)%cycles(2) - S(j)%cycles(1)) / S(j)%cycles(3)
		i = i + S(j)%cycles(2) - S(j)%cycles(1) + 1
	else
		write (*,610) S(j)%sat,S(j)%cycles(1:2),pass
	endif
enddo
if (ptrx1 > ntrx) then
	write (*,615) 'average of all cycles',i
	i = i + 1
endif
if (ptrx1 > ntrx+1) write (*,615) 'standard deviation of all cycles',i

! Describe variables
write (*,620)
i = 1
do j = 1,nsel
	write (*,622) i,i+ptrx1-ptrx0,trim(S(nsat)%sel(j)%info%long_name),trim(S(nsat)%sel(j)%info%units)
	i = i + ptrx1 - ptrx0 + 1
enddo
if (counter) then
	write (*,625) i,'number of measurements'
	i = i + 1
endif
if (numbered) write (*,625) i,'record number'
write (*,630)

end subroutine write_header

!***********************************************************************

end program radscolin
