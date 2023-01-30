      subroutine cmsave
      integer k,l
      include 'xrgb.inc'

      open (20,file=cmfile)
      rewind (20)
      write (20,550) '#COL'
      write (20,550) '#'
      write (20,550) '# This file was created by xrgb'
      write (20,550) '#'
      if (cmap(1,0).ge.0) write (20,550) 'C_BG   ',(cmap(l,0),l=1,3)
      if (cmap(1,1).ge.0) write (20,550) 'C_FG   ',(cmap(l,1),l=1,3)
      if (cmap(1,2).ge.0) write (20,550) 'C_LAND ',(cmap(l,2),l=1,3)
      if (cmap(1,3).ge.0) write (20,550) 'C_BAD  ',(cmap(l,3),l=1,3)
      if (cmap(1,4).ge.0) write (20,550) 'C_COAST',(cmap(l,4),l=1,3)
      if (cmap(1,5).ge.0) write (20,550) 'C_CONT ',(cmap(l,5),l=1,3)
      l=0
      do k=6,maxind
	 if (cmap(1,k).ge.0) l=l+1
      enddo
      write (20,550) '#'
      write (20,550) '# Now follows a colour range of',l
      write (20,550) '#'
      do k=6,maxind
	 if (cmap(1,k).ge.0) write (20,550) 'C_RNG   ',(cmap(l,k),l=1,3)
      enddo
      close (20)
      changed=.false.
550   format (a,3i4)
      end
