!-----------------------------------------------------------------------
! $Id$
!
! Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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

contains

!*erspass - Determine orbit nr, phase, cycle nr and pass nr for ERS-1/2
!+
logical function erspass (ers, utc, orbitnr, phasenm, cyclenr, passnr, tnode, lnode)
use rads
real(eightbytereal), intent(in) :: utc
integer(fourbyteint), intent(in) :: ers
integer(fourbytereal), intent(out) :: orbitnr, cyclenr, passnr
real(eightbytereal), intent(out) :: tnode, lnode
character(1), intent(out) :: phasenm

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
integer(fourbyteint), parameter :: mpass=100000
type :: passtable
	integer(fourbyteint) :: orbitnr
	character(1) :: phasenm
	integer(fourbyteint) :: cyclenr,passnr
	real(eightbytereal) :: start,tnode,lnode
endtype
type(passtable) :: q(mpass)
character(160) :: line

save pnt,npass,olders,q

! Load the file with pass definitions upon initialisation

if (ers == olders) then
	new = .false.
else
	write (line,610) trim(radsdataroot),ers,ers
	unit = freeunit()
	npass=0
	open (unit,file=line,status='old')
	do
		read (unit,'(a)',iostat=ios) line
		if (ios /= 0) exit
		if (line(:1) == '#') cycle
		npass = npass + 1
		if (npass > mpass) call fin('erspass: too many passes')
		read (line,600) q(npass)
	enddo
	close (unit)
	olders = ers
	new = .true.
	pnt = 1
endif
600 format (i6,1x,a1,i4,i5,2f14.3,f11.6)
610 format (a,'/../ext/e',i1,'/ers',i1,'passes.tab')

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

end module rads_devel
