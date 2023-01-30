      subroutine noneg(kb)
      implicit none
      integer kb
*
* Print message and convert kb to positive
*
  550 format (a)
      if (kb.lt.0) then
        write (0,550) 'ingrid: negative not allowed, abs() used'
        kb=abs(kb)
      endif
      end
