      subroutine cmcopy
      character*80 text
      include 'xrgb.inc'
      integer idx(2),m

1     call menu(0,' ')
      call menu(12,'Done')
      call menu(13,'Click colour to copy')
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
      write (text,700) 'Copy colour index ',ctxt(idx(1)),' to ...'
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
      write (text,700) 'Copying colour index ',ctxt(idx(1)),
     .' to ',ctxt(idx(2))
      call menu(13,text)

* Do the actual copying

      changed=.true.
      do m=1,3
         cmap(m,idx(2))=cmap(m,idx(1))
      enddo

* Update new colour

      call colupdt(idx(2))

* Remove the tag

      call untag(idx(1))
      call untag(idx(2))
      goto 1
      end
