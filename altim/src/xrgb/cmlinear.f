      subroutine cmlinear
      character*80 text
      integer idx(2)
      include 'xrgb.inc'
      real dx,dy,a,b
      integer k,l

1     call menu(0,' ')
      call menu(12,'Done')
      call menu(13,
     .'Linear interpolation of colours. Click first colour.')
2     call pgcurse(x,y,ch)
      call whatindx(idx(1))
      ich=ichar(ch)
      if (ich.ne.65) then
	 goto 2
      else if (idx(1).lt.0) then
	 return
      endif
      call tag(idx(1))
700   format(a,a,a,a)
      write (text,700) 'Linear interpolation from ',ctxt(idx(1)),
     .' to ...'
      call menu(13,text)

3     call pgcurse(x,y,ch)
      call whatindx(idx(2))
      ich=ichar(ch)
      if (ich.ne.65) then
	 goto 3
      else if (idx(2).lt.0) then
	 call untag(idx(1))
         return
      endif
      call tag(idx(2))
      write (text,700) 'Linear interpolation from ',ctxt(idx(1)),
     .' to ',ctxt(idx(2))
      call menu(13,text)
c
c Do the actual linear thing
c
      changed=.true.
      dx=idx(2)-idx(1)
      do k=idx(1),idx(2)
         do l=1,3
c
c Determine tilt and bias
c
            dy=cmap(l,idx(2))-cmap(l,idx(1))
            a=dy/dx
            b=cmap(l,idx(1))-a*idx(1)
            cmap(l,k)=nint(a*k+b)
	 enddo
	 call colupdt(k)
      enddo

      call untag(idx(1))
      call untag(idx(2))
      goto 1
      end
