**PIMOPT -- Strip options from the command cards
*+
      FUNCTION PIMOPT (FLAG, CMD, VAL1, VAL2, VAL3, VAL4)
      LOGICAL PIMOPT
      CHARACTER*(*) FLAG, CMD
      REAL VAL1, VAL2, VAL3, VAL4
*
* This routine scans the command cards for options starting with FLAG.
* Then it reads the values after the flag and strips of the option from
* the line. It returns PIMOPT=.TRUE. when the option is used.
*
      integer lf,lo,i,ios

      pimopt=.false.

      i=index(cmd,flag)
      if (i.le.0) return

      lf=len(flag)
      lo=index(cmd(i:),' ')
      if (lf.le.lo-2) then
         read (cmd(i+lf:i+lo-2),*,iostat=ios) val1,val2,val3,val4
         cmd(i:)=cmd(i+lo:)
      endif
      pimopt=.true.
      end

**DPIMOPT -- Strip options from the command cards
*+
      FUNCTION DPIMOPT (FLAG, CMD, VAL1, VAL2, VAL3, VAL4)
      LOGICAL DPIMOPT
      CHARACTER*(*) FLAG, CMD
      REAL*8 VAL1, VAL2, VAL3, VAL4
*
* This routine scans the command cards for options starting with FLAG.
* Then it reads the values after the flag and strips of the option from
* the line. It returns DPIMOPT=.TRUE. when the option is used.
*
      integer lf,lo,i,ios

      dpimopt=.false.

      i=index(cmd,flag)
      if (i.le.0) return

      lf=len(flag)
      lo=index(cmd(i:),' ')
      if (lf.le.lo-2) then
         read (cmd(i+lf:i+lo-2),*,iostat=ios) val1,val2,val3,val4
         cmd(i:)=cmd(i+lo:)
      endif
      dpimopt=.true.
      end
