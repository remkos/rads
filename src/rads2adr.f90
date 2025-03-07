!-----------------------------------------------------------------------
! Copyright (c) 2011-2025  Remko Scharroo
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

!*rads2adr -- Select RADS data and output to ADR or XGF
!+
program rads2adr
!
! This program converts the RADS NetCDF altimeter data to the old
! DEOS formats aADR, xADR and XGF, depending on whether the program
! is called as rads2adr, rads2xadr or rads2xgf.
! At the same time it applies the standard selection criteria
! and allows some further modifications and selections.
!-----------------------------------------------------------------------
use rads
use rads_netcdf
use rads_misc
use rads_time

! RADS structures
type(rads_sat) :: S
type(rads_pass) :: P

! Local declarations, etc.
integer(fourbyteint) :: outunit
character(len=rads_cmdl) :: outname = '', dirname = ''
integer(fourbyteint) :: i, l, ios, cycle, pass, step = 1, nseltot = 0, nselmax = huge(0_fourbyteint)
integer(fourbyteint) :: reject = -1

! Operation modes
character(len=80) :: prog
integer, parameter :: mode_adr = 2, mode_xadr = 3, mode_xgf = 4
integer :: mode

! Output file definitions
integer(fourbyteint) :: adr4(5), xgf4(4), rec_length, nselfile = 0, nselpass = 0
integer(twobyteint) :: adr2(4), xgf2
character(len=4) :: suffix

! Initialize RADS or issue help
call synopsis
call rads_set_options ('o:r:: output: reject-on-nan:: maxrec: step:')
call rads_init (S)
if (S%error /= rads_noerr) call rads_exit ('Fatal error')

! Set default options related to output data type
do i = 1,rads_nopt
	select case (rads_opt(i)%opt)
	case ('o', 'output')
		l = len_trim(rads_opt(i)%arg)
		if (rads_opt(i)%arg(l:l) == '/') then
			dirname = rads_opt(i)%arg
		else
			outname = rads_opt(i)%arg
		endif
	case ('r', 'reject-on-nan')
		call rads_parse_r_option (S, rads_opt(i)%opt, rads_opt(i)%arg, reject)
	case ('maxrec')
		read (rads_opt(i)%arg, *, iostat=ios) nselmax
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	case ('step')
		read (rads_opt(i)%arg, *, iostat=ios) step
		if (ios /= 0) call rads_opt_error (rads_opt(i)%opt, rads_opt(i)%arg)
	end select
enddo

! Force the -V option to conform to the operation mode (rads2adr, rads2xadr, rads2xgf)
if (S%nsel > 0) then
	call rads_error (S, rads_noerr, 'Ignoring -V|--var|--sel option')
	S%nsel = 0
	deallocate (S%sel)
endif

! Determine operation mode
call getarg (0, prog)
if (index(prog, 'rads2adr') > 0) then
	mode = mode_adr
	call rads_parse_varlist (S, '1,2,3,4,0,20')
	suffix = '.adr'
	rec_length = 24
else if (index(prog, 'rads2xadr') > 0) then
	mode = mode_xadr
	call rads_parse_varlist (S, '1,2,3,4,0,17,10,5')
	suffix = '.adr'
	rec_length = 28
else
	mode = mode_xgf
	call rads_parse_varlist (S, '1,2,3,0,20')
	suffix = '.xgf'
	rec_length = 18
endif

! If SLA is among the results, remember which index that is
do i = 1,S%nsel
	if (S%sel(i)%info%datatype == rads_type_sla) then
		if (reject == -1) reject = i
	endif
enddo

! Open single output file, if requested
if (outname /= '') then
	outunit = getlun()
	open (outunit, file=outname, access='direct', form='unformatted', status='replace', recl=rec_length)
endif

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

! Close data file before exit
if (outname /= '') call close_datafile

! Print overall statistics and close RADS
call rads_stat (S)
call rads_end (S)

contains

!***********************************************************************

subroutine synopsis
if (rads_version ('Select RADS altimeter data and output to ADR or XGF')) return
call rads_synopsis
write (*,1300)
1300 format (/ &
'Program specific [program_options] are:'/ &
'  -r, --reject-on-nan VAR   Reject records if variable VAR on -V specifier is NaN'/ &
'  -r #                      Reject lines if data item number # on -V specifier is NaN'/ &
'  -r 0, -r none, -r         Do not reject lines with NaN values'/ &
'  -r n, -r any              Reject lines if any value is NaN'/ &
'                      Note: If no -r option is given -r sla is assumed'/ &
'  --step N                  Step through records with stride N (default = 1)'/ &
'  -o, --output OUTNAME      Specify name of a single output file or (when ending in /) directory'/ &
'                            name for pass files; default is pass files in current directory')
stop
end subroutine synopsis

!***********************************************************************
! Process data for a single pass

subroutine process_pass (ndata, nsel)
integer(fourbyteint), intent(in) :: ndata, nsel
real(eightbytereal) :: data(ndata,nsel)
character(len=80) :: passname
integer(fourbyteint) :: i

! Open a new output file when outname is not specified
if (outname == '') then
	write (passname, '(a,a2,"p",i4.4,"c",i3.3,a)') trim(dirname), S%sat, cycle, pass, suffix
	outunit = getlun()
	open (outunit, file=passname, access='direct', form='unformatted', status='replace', recl=rec_length)
	nselfile = 0
endif
nselpass = 0

! Read the data
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
	nselfile = nselfile + 1

	! Write out the data in the appropriate format
	if (mode == mode_adr) then
		adr4(1) = floor(data(i,1))
		adr4(2) = nint((data(i,1)-adr4(1))*1d6)
		adr4(3:4) = nint4(data(i,2:3)*1d6)
		adr4(5) = nint4(data(i,4)*1d3-700d0)
		adr2(1:2) = nint2(data(i,5:6)*1d3)
		write (outunit,rec=nselfile+1) adr4,adr2(1:2)
	else if (mode == mode_xadr) then
		adr4(1) = floor(data(i,1))
		adr4(2) = nint((data(i,1)-adr4(1))*1d6)
		adr4(3:4) = nint4(data(i,2:3)*1d6)
		adr4(5) = nint4(data(i,4)*1d3-700d0)
		adr2(1:4) = nint2(data(i,5:8)*1d3)
		write (outunit,rec=nselfile+1) adr4,adr2(1:4)
	else
		xgf4(1) = nint(data(i,1))
		xgf4(2:4) = nint4(data(i,2:4)*1d6)
		xgf2 = nint2(data(i,5)*1d3)
		write (outunit,rec=nselfile+1) xgf4,xgf2
	endif
enddo

nseltot = nseltot + nselpass

! Close per-pass output file
if (outname /= '') then
	! Keep file open
else if (nselfile == 0) then
	close (outunit,status='delete')
else
	call close_datafile
endif

end subroutine process_pass

!***********************************************************************

subroutine close_datafile
integer(twobyteint), parameter :: extra(2) = (/3_twobyteint,0_twobyteint/)
integer(twobyteint) :: bound(4)
bound(1:2) = nint2(S%lon%info%limits)
bound(3:4) = nint2(S%lat%info%limits)
if (mode == mode_adr) then
	write (outunit,rec=1) 'aADR', S%satellite(1:8), bound, nselfile
else if (mode == mode_xadr) then
	write (outunit,rec=1) 'xADR', S%satellite(1:8), bound, nselfile, extra
else
	write (outunit,rec=1) '@XGF', nselfile
endif
close (outunit, status='keep')
end subroutine close_datafile

!***********************************************************************

end program rads2adr
