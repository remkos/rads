      character*80 c
      read (5,'(a)') c
      read (c,*,iostat=ios) i,j
      write (6,*) i,j,ios
      end

