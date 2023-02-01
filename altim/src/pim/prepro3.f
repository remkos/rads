      subroutine prepro3(file)

* Third level preprocessing of setup file.
* To avoid recursive programming we have prepro1, prepro2, prepro3,
* with the same use.

      character argum*1024,filenm*1024,command*5,file*(*)
      logical exist
      integer l,lnblnk,unit,freeunit

      if (file.eq.'-') then
	 unit=5
	 write (*,610) '(stdin)'
      else
	 inquire (file=file,exist=exist)
	 if (exist) then
	    filenm=file
	 else
	    filenm='/user/altim'
	    call checkenv('ALTIM',filenm,l)
	    filenm(l+1:)='/pim/'//file
	 endif
	 unit=freeunit()
	 write (*,610) filenm(:lnblnk(filenm))
	 open (unit,file=filenm,status='old')
      endif
10    read (unit,550,end=999) argum
      if (argum.eq.' ') goto 10
      if (argum(:1).eq.'#') goto 10
      call strip(argum,command)
      call grtoup(command,command)
      if (argum(:3).eq.'REM') goto 10
      l=lnblnk(argum)
      if (command.eq.'INCLU' .or. command.eq.'COLOU') then
	 stop 'Too many include levels'
      else
	 write (7,600) command,argum(:l)
      endif
      goto 10

550   format (a)
600   format (a5,1x,a)
610   format ('... reading ',a)

999   if (unit.ne.5) call grflun(unit)
      end
