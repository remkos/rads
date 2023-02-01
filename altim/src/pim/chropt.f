**CHROPT -- Strip options from the command cards
*+
      FUNCTION CHROPT (FLAG, CMD, STRING)
      LOGICAL CHROPT
      CHARACTER*(*) FLAG, CMD, STRING
*
* This routine scans the command cards for options starting with FLAG.
* Then it reads the string after the flag and strips of the option from
* the line. It returns CHROPT=.TRUE. when the option is used.
*
      integer lf,lo,i

      chropt=.false.

      i=index(cmd,flag)
      if (i.le.0) return

      lf=len(flag)
      lo=index(cmd(i:),' ')
      string=cmd(i+lf:i+lo-2)
      cmd(i:)=cmd(i+lo:)
      chropt=.true.
      end
