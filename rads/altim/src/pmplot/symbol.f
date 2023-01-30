      program symbol

      real vpdef(4),x,y
      integer k,l,m
      character mchar*4,text*80
      data vpdef /.05,.95,.05,.95/

      call getarg(1,text)
      if (text.eq.' ') text='?'
      call pgbeg(0,text,1,1)
      call pgsvp(vpdef(1),vpdef(2),vpdef(3),vpdef(4))
      call pgswin(0.,15.,0.,10.)
      call pgbox('BCG',1.,0,'BCG',1.,0)
   50 format(i4.4)
      call getarg(2,text)
      if (text.eq.' ') text='750'
      read (text,*) m
      m=m-1
      do l=10,1,-1
         do k=1,15
            m=m+1
            write (mchar,50) m
            x=k-.5
            y=l-.5
            call pgsch(2.)
            call pgpt(1,x,y,m)
            call pgsch(.5)
            x=k-.95
            y=l-.13
            call pgptxt(x,y,0.,0.,mchar)
         enddo
      enddo
      call pgend
      end
