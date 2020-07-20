!-----------------------------------------------------------------------
! Copyright (c) 2011-2020  Remko Scharroo
! See LICENSE.TXT file for copying and redistribution conditions.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!-----------------------------------------------------------------------

!*rads2nc -- Select RADS data and output to NetCDF
!+
program rads2nc
!
! This program converts the RADS NetCDF altimeter data and output them
! to new NetCDF data files.
! At the same time it applies the standard selection criteria
! and allows some further modifications and selections.
!-----------------------------------------------------------------------
use rads
use rads_netcdf
use rads_misc
use netcdf

! RADS structures
type(rads_sat) :: S
type(rads_pass) :: P, Pout

! Local declarations, etc.
character(len=rads_cmdl) :: outname = '', dirname = ''
integer(fourbyteint) :: i, l, ios, cycle, pass, step = 1, nseltot = 0, nselmax = huge(0_fourbyteint)
integer(fourbyteint) :: reject = -1

! Output file definitions
integer(fourbyteint) :: nselpass = 0

! Initialize RADS or issue help
call synopsis
call rads_set_options ('r::o: reject-on-nan:: output: step: maxrec:')
call rads_init (S)
if (S%error /= rads_noerr) call rads_exit ('Fatal error')

! Scan command line arguments
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('o', 'output')
		l = len_trim(rads_opt(i)%arg)
		if (rads_opt(i)%arg(l:l) == '/') then
			dirname = rads_opt(i)%arg
		else
			outname = rads_opt(i)%arg
		endif
	case ('r', 'recect-on-nan')
		call rads_parse_r_option (S, rads_opt(i)%opt, rads_opt(i)%arg, reject)
	case ('maxrec')
		read (rads_opt(i)%arg, *, iostat=ios) nselmax
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('step')
		read (rads_opt(i)%arg, *, iostat=ios) step
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	end select
enddo

! Check if there is at least one variable specified
if (S%nsel == 0) call rads_exit ('No variables are selected')

! If no -r option was used and SLA is among the selected variables, remember its index
do i = 1,S%nsel
	if (reject == -1 .and. S%sel(i)%info%datatype == rads_type_sla) reject = i
enddo

! Initialise an output pass structure
call rads_init_pass_struct (S, Pout)

! Now loop through all cycles and passes
do cycle = S%cycles(1), S%cycles(2), S%cycles(3)
	! Stop processing after too many output lines
	if (nseltot >= nselmax) then
		write (*,760) nseltot,nselmax
		exit
	endif

	! Process passes one-by-one
	do pass = S%passes(1), S%passes(2), S%passes(3)
		call rads_open_pass (S, P, cycle, pass)
		if (P%ndata > 0) call process_pass (P%ndata, S%nsel)
		if (rads_verbose >= 1) call rads_progress_bar (S, P, nselpass)
		call rads_close_pass (S, P)
	enddo
enddo
760 format(/'Maximum number of output records reached (',i0,' >= ',i0,')')

! Finish progress bar
if (rads_verbose >= 1) write (*,*)

! Close data file before exit and/or empty out history
call rads_close_pass (S, Pout)

! Print overall statistics and close RADS
call rads_stat (S)
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('Select RADS altimeter data and output to NetCDF')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r, --reject-on-nan VAR   Reject records if variable VAR on -V specifier is NaN'/ &
'  -r NR                     Reject records if data item number NR on -V specifier is NaN'/ &
'  -r 0, -r none, -r         Do not reject records with NaN values'/ &
'  -r n, -r any              Reject records if any value is NaN'/ &
'                      Note: If no -r option is given -r sla is assumed'/ &
'  --step N                  Step through records with stride n (default: 1)'/ &
'  --maxrec NREC             Specify maximum number of output records (default: unlimited)'/ &
'  -o, --output OUTNAME      Specify name of a single output file or (when ending in /) directory'/ &
'                            name for pass files; default is pass files in current directory')
stop
end subroutine synopsis

!***********************************************************************

subroutine process_pass (ndata, nsel)
integer(fourbyteint), intent(in) :: ndata, nsel
real(eightbytereal) :: data(ndata,nsel)
logical :: accept(ndata)
integer(fourbyteint) :: i, start(1)

! Read the data
nselpass = 0
accept = .false.
do i = 1,nsel
	call rads_get_var (S, P, S%sel(i), data(:,i))
enddo

! Loop through the data
do i = 1,ndata,step
	! See if we have to reject this record
	if (reject > 0) then
		if (isnan_(data(i,reject))) cycle
	else if (reject == -2) then
		if (any(isnan_(data(i,:)))) cycle
	endif
	nselpass = nselpass + 1
	accept(i) = .true.
enddo

! If no data left, then return
if (nselpass == 0) return

! Open output pass file
start(1) = nseltot + 1
if (outname == '') then
	Pout = P
	nullify (Pout%history)
	call rads_create_pass (S, Pout, nselpass, name=dirname)
	call rads_def_var (S, Pout, S%sel)
	start(1) = 1
else if (nseltot == 0) then
	call rads_create_pass (S, Pout, 0, name=outname)
	call rads_def_var (S, Pout, S%sel)
endif

! Write the variables
do i = 1,nsel
	call rads_put_var (S, Pout, S%sel(i), pack(data(:,i),accept), start)
enddo
nseltot = nseltot + nselpass

! Close per-pass output file
! We need to keep the history, etc, since rads_close_pass (S, P) in the main program
! is going to deallocate those
if (outname == '') call rads_close_pass (S, Pout, .true.)

end subroutine process_pass

!***********************************************************************

end program rads2nc
