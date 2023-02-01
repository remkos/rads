**FIN -- Routine to end or abort program gracefully
*+
      SUBROUTINE FIN (STRING)
      CHARACTER*(*) STRING

* This routine terminates the execution of a program gracefully.
* If STRING is empty, the program ends normally with the message
* "Normal end of program xxx" printed to standard output.
* Otherwise it first echos STRING and then "Execution of program
* xxx terminated".
* After that the STOP statement is excuted in both cases.
*
* Argument:
*  STRING (input): Message echoed at program termination. If empty, the
*                  program ends normally.
*-
* 10-Aug-1998 -- Created by Remko Scharroo
*-----------------------------------------------------------------------
      character*80 prog
      integer l,lnblnk

      if (string.eq.' ') then
         call getarg(0,prog)
         l=lnblnk(prog)
         write (*,550) 'Normal end of program '//prog(:l)
      else
         l=lnblnk(string)
         write (*,550) string(:l)
         call getarg(0,prog)
         l=lnblnk(prog)
         write (*,550) 'Execution of program '//prog(:l)//
     |		' terminated'
      endif
      stop
550   format(a)
      end
