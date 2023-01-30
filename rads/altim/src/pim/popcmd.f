**POPCMD -- Pop last matching command from strack
*
      function popcmd(cmd,argum)
      logical popcmd
      character*(*) cmd,argum
      character command*5,argument*1024

      popcmd=.false.
      rewind(7)

10    read (7,600,end=90) command,argument
      if (command.eq.cmd) then
         argum=argument
         popcmd=.true.
      endif
      goto 10

90    rewind(7)
600   format (a5,1x,a)
      end

**POP1CMD -- Pop first matching command from strack
*
      function pop1cmd(cmd,argum)
      logical pop1cmd
      character*(*) cmd,argum
      character command*5,argument*1024

10    read (7,600,end=90) command,argument
      if (command.ne.cmd) goto 10
      argum=argument
      pop1cmd=.true.
      return

90    pop1cmd=.false.
      rewind(7)
600   format (a5,1x,a)
      end

**POPNEXTCMD -- Pop first matching command from strack
*
      function popnextcmd(cmd,argum)
      logical popnextcmd
      character*(*) cmd,argum
      character command*5,argument*1024

      read (7,600,end=90) command,argument
      backspace (7)
      if (command.ne.cmd) goto 90
      argum=argument
      popnextcmd=.true.
      return

90    popnextcmd=.false.
600   format (a5,1x,a)
      end
