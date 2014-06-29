!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2014  Remko Scharroo (Altimetrics LLC)
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
use typesizes
use rads_time
use rads_misc

contains

!*erspass - Determine orbit nr, phase, cycle nr and pass nr for ERS-1/2
!+
logical function erspass (ers, utc, orbitnr, phasenm, cyclenr, passnr, tnode, lnode)
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
integer(fourbyteint) :: unit,freeunit,npass=0,pnt,olders=0,ios
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
	unit = freeunit()
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
else
	write (*,1301) trim(arg),syntax
endif
1300 format (/a,' [rads_dataselectors]',a // &
'Optional [rads_dataselectors] are:' / &
'  -C, --cycle=C0[,C1]       Select data for one or more cycles' / &
'  --time=T0,T1              Specify time selection (optionally use --ymd=, --doy=,' / &
'                            or --sec= for [YY]YYMMDD[HHMMSS], YYDDD, or SEC85)')
1301 format (/a,' [rads_dataselectors]',a // &
'Required argument is:' / &
'  -S, --sat=SAT[/PHASE]     Specify satellite [and phase] (e.g. e1/g, tx)' // &
'Optional [rads_dataselectors] are:' / &
'  -C, --cycle=C0[,C1[,DC]]  Specify first and last cycle and modulo' / &
'  -P, --pass=P0[,P1[,DP]]   Specify first and last pass and modulo' / &
'  --time=T0,T1              Specify time selection (optionally use --ymd=, --doy=,' / &
'                            or --sec= for [YY]YYMMDD[HHMMSS], YYDDD, or SEC85)')
end subroutine synopsis_devel
end module rads_devel
