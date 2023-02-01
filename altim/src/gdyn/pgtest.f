      character*4 ch
      real x,y

      call pgbeg(0,"/xs",1,1)
      call pgswin(0.0,1.0,0.0,1.0)
      call pgbox('abcnst',0.0,0,'abcnst',0.0,0)

10    call pgcurs(x,y,ch)
      write (6,'(2f6.3,i5)') x,y,ichar(ch)
      goto 10
      end
