!-----------------------------------------------------------------------
! Copyright (c) 2011-2018  Remko Scharroo
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

module rads_devel

contains

!*erspass - Determine orbit nr, phase, cycle nr and pass nr for ERS-1/2
!+
logical function erspass (ers, utc, orbitnr, phasenm, cyclenr, passnr, tnode, lnode)
use rads_misc
real(eightbytereal), intent(in) :: utc
integer(fourbyteint), intent(in) :: ers
integer(fourbytereal), intent(out) :: orbitnr, cyclenr, passnr
real(eightbytereal), intent(out) :: tnode, lnode
character(len=1), intent(out) :: phasenm

! This routine returns the ESA orbit number (orbitnr), mission phase
! (phase), and the cycle number (cyclenr) and pass (or half-orbit) number
! (passnr) for any given time during the ERS-1 or ERS-2 mission.
!
! Upon initialisation one of the files $(RADSROOT)/ext/e?/ers?passes.tab
! is loaded.
!
! Input arguments:
!  ers     : 1 = ERS-1, 2 = ERS-2
!  utc     : Time in UTC seconds from 1985 (SEC85)
!
! Output arguments:
!  orbitnr : ESA orbit number
!  phasenm : Mission phase (A-G)
!  cyclenr : Cycle number (1-156)
!  passnr  : Pass number (1-4882)
!  tnode   : Time of the equator crossing (SEC85)
!  lnode   : Longitude of the equator crossing (deg)
!  erspass : .TRUE. if the pass has changed
!-----------------------------------------------------------------------
logical :: new
integer(fourbyteint) :: unit,npass=0,pnt,olders=0,ios
integer(fourbyteint), parameter :: mpass=170000, mjd90=1826 ! Days from 1985 to 1990
real(eightbytereal) :: mjd
type :: passtable
	integer(fourbyteint) :: orbitnr
	character(len=1) :: phasenm
	integer(fourbyteint) :: cyclenr,passnr
	real(eightbytereal) :: start,tnode,lnode
endtype
type(passtable) :: q(mpass)
character(len=160) :: line

save pnt,npass,olders,q

! Load the file with pass definitions upon initialisation

if (ers == olders) then
	new = .false.
else
	call parseenv ('${RADSROOT}/ext/reaper/', line)
	write (line,610) trim(line), ers
	unit = getlun()
	npass = 0
	open (unit, file=line, status='old')
	read (unit, *, iostat=ios) ! Skip header
	do
		read (unit,'(a)',iostat=ios) line
		if (ios /= 0) exit

		npass = npass + 1
		read (line,600) mjd
		q(npass)%start = (mjd + mjd90) * 86400d0

		read (unit,'(a)',iostat=ios) line
		if (ios /= 0) exit

		read (line,600) mjd, q(npass)%phasenm, q(npass)%cyclenr, q(npass)%passnr, q(npass)%orbitnr, q(npass)%lnode
		q(npass)%tnode = (mjd + mjd90) * 86400d0
	enddo
	close (unit)
	olders = ers
	new = .true.
	pnt = 1
endif
600 format (f11.6,21x,a1,i4,i5,i6,2f8.3)
610 format (a,'ER',i1,'_ORF.txt')

! Look for the table entry based on the utc time

do
	if (utc > q(pnt+1)%start .and. pnt < npass) then
		new = .true.
		pnt = pnt + 1
	else if (utc < q(pnt)%start .and. pnt > 1) then
		new = .true.
		pnt = pnt - 1
	else
		exit
	endif
enddo

! Return pass information

phasenm = q(pnt)%phasenm
orbitnr = q(pnt)%orbitnr
cyclenr = q(pnt)%cyclenr
passnr = q(pnt)%passnr
tnode = q(pnt)%tnode
lnode = q(pnt)%lnode

erspass = new

end function erspass

!*synopsis_devel -- Provide synopsis shared by development tools
!+
subroutine synopsis_devel (syntax)
character(len=*), intent(in) :: syntax
!
! This routine produces the usage statement plus some generally
! used options for the development tools. A distinction is made between
! rads_gen_* commands and others.
!
! Argument:
!   syntax : string to be added to syntax
!-----------------------------------------------------------------------
character(len=160) :: arg
call getarg (0, arg)
if (index(arg, 'rads_gen_') > 0) then
	write (*,1300) trim(arg),syntax
	write (*,1310)
else
	write (*,1301) trim(arg),syntax
	write (*,1310)
endif
1300 format (/'Usage: ',a,' [rads_dataselectors]',a // &
'Optional [rads_dataselectors] are:' / &
'  -S, --sat SAT[/PHASE]     Specify satellite [and phase] (e.g. e1/g, tx)' / &
'  -C, --cycle C0[,C1]       Select data for one or more cycles')
1301 format (/'Usage: ',a,' [rads_dataselectors]',a // &
'Required argument is:' / &
'  -S, --sat SAT[/PHASE]     Specify satellite [and phase] (e.g. e1/g, tx)' // &
'Optional [rads_dataselectors] are:' / &
'  -C, --cycle C0[,C1[,DC]]  Specify first and last cycle and modulo' / &
'  -P, --pass P0[,P1[,DP]]   Specify first and last pass and modulo')
1310 format ( &
'  --time T0,T1              Specify time selection (optionally use --ymd, --doy,' / &
'                            or --sec for [YY]YYMMDD[HHMMSS], YYDDD, or SEC85)' // &
'Common [rads_options] are:'/ &
'  --help                    Print this syntax massage'/ &
'  --log FILENAME            Send statistics to FILENAME (default is standard output)'/ &
'  -q, --quiet               Suppress warning messages (but keeps fatal error messages)' / &
'  -v, --verbose             Increase verbosity level'/ &
'  --debug LEVEL             Set debug/verbosity level'/ &
'  --version                 Version info')
end subroutine synopsis_devel

!*log_string -- Print string to log output
!+
subroutine log_string (string, advance)
use rads
character(len=*), intent(in) :: string
logical, intent(in), optional :: advance
!
! This routine prints to the log output (rads_log_unit) a trimmed string
! followed by ' ... ' or by a carriage return if <advance> is true.
!
! Argument:
!   string : string to be printed
!   advance: (optional) add carriage return if .true. (default is .false.)
!-----------------------------------------------------------------------
550 format (a)
552 format (a, ' ... ')
if (present(advance) .and. advance) then
	write (rads_log_unit,550) trim(string)
else
	write (rads_log_unit,552,advance='no') trim(string)
endif
end subroutine log_string

!*log_pass -- Print filename of pass to log output
!+
subroutine log_pass (P, advance)
use rads
type(rads_pass), intent(in) :: P
logical, intent(in), optional :: advance
!
! This routine prints to the log output (rads_log_unit) the name of the
! pass file, followed by ' ... ' by a carriage return if <advance> is true.
!
! Argument:
!   P      : pass structure
!   advance: (optional) add carriage return if .true. (default is .false.)
!-----------------------------------------------------------------------
call log_string (P%fileinfo(1)%name(len_trim(P%S%dataroot)+2:), advance)
end subroutine log_pass

!*log_records -- Print string to log output
!+
subroutine log_records (count, P)
use rads
integer(fourbyteint), intent(in) :: count
type(rads_pass), intent(in), optional :: P
!
! This routine prints to the log output (rads_log_unit) one of two
! strings, depending on whether <P> is supplied.
! Without <P>: xxxx records changed
! With <P>   : xxxx records written to <pass_file_name>
! The <pass_file_name> is trimmed and the output is followed by a
! carriage return.
!
! Argument:
!   count  : number to be printed
!   string : string to be printed (optional)
!-----------------------------------------------------------------------
551 format (i4,' records written to ',a)
553 format (i4,' records changed')

if (present(P)) then
	write (rads_log_unit,551) count, trim(P%fileinfo(1)%name(len_trim(P%S%dataroot)+2:))
else
	write (rads_log_unit,553) count
endif
end subroutine log_records

end module rads_devel
