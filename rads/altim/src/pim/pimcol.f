**PIMCOL -- Read and set colour range (with or without shading)
*+
      SUBROUTINE PIMCOL (SHADE)
      LOGICAL SHADE
*
* This routine reads and sets the colour representation for the colour range.
* If requested, it shades the colours.
*
* Important variables
* NC1 = Number of colors in colour range - 1
* NC2 = Number of shades for each colour - 1
*-
      include 'pim.inc'

      integer nquant,i,j,ic
      real h,l,s,r,g,b

* Get device colour capabilities

      call pgqcol (i,c_1)

* Read and set specified colours

      if (popcmd('C_BG ',argum)) call setcol(c_bg   ,argum)
      if (popcmd('C_FG ',argum)) call setcol(c_fg   ,argum)
      if (popcmd('C_LAN',argum)) call setcol(c_land ,argum)
      if (popcmd('C_BAD',argum)) call setcol(c_bad  ,argum)
      if (popcmd('C_COA',argum)) call setcol(c_coast,argum)
      if (popcmd('C_CON',argum)) call setcol(c_cont ,argum)
      if (popcmd('C_BLU',argum)) call setcol(c_blue ,argum)
      if (popcmd('C_RED',argum)) call setcol(c_red  ,argum)
      if (popcmd('C_GRE',argum)) call setcol(c_grey ,argum)
      if (popcmd('C_LGR',argum)) call setcol(c_lgrey,argum)
      if (popcmd('C_DGR',argum)) call setcol(c_dgrey,argum)
      if (popcmd('C_WHI',argum)) call setcol(c_white,argum)
      if (popcmd('C_BLA',argum)) call setcol(c_black,argum)
      rewind(7)

      if (postscript) then
         call setcol(c_bg,'255 255 255')
         call setcol(c_fg,'000 000 000')
      endif

* Read and set colour range

      ic=c_0
10    if (pop1cmd('C_RNG',argum)) then
	 if (ic.gt.c_1)
     .      stop 'Number of colours exceeds device capabilities'
	 call setcol(ic,argum)
	 ic=ic+1
         goto 10
      endif
      nc1=(ic-1)-c_0

      if (shade) then

* Use maximum number of indices as default for this device. For true
* color devices (i<0) use nc2=255

	 nquant=c_1+1
	 nc2=255
	 if (i.lt.0) nquant=(nc1+1)*(nc2+1)+c_0
	 if (popcmd('QUANT',argum)) then
	    nquant=c_1+1
	    read (argum,*,iostat=i) nquant
	 endif
	 nc2=(nquant-c_0)/(nc1+1)-1
	 write (0,600) nc1+1,nc2+1
	 if (nc2.gt.255)
     .      stop 'Number of requested shades exceeds capabilities'
	 ic=c_0

* Set shades for each colour. Do not use black and white, but use
* lightnesses between 0.20 and 0.95.
* Set color indices up to C_1,  then use HEXcodes (true color devices only).

	 do i=0,nc1
	    call grxhls(rgb(1,c_0+i),rgb(2,c_0+i),rgb(3,c_0+i),h,l,s)
 	    do j=0,nc2
	       l=shade0+(shade1-shade0)*j/(nc2+1)
	       call grxrgb(h,l,s,r,g,b)
*	       if (ic.le.c_1) then
		  call pgscr(ic,r,g,b)
		  sh_mat(c_0+i,j)=ic
* GRXHEX not used anymore. NO HEXcodes
*	       else
*		  call grxhex(r,g,b,ic)
*		  sh_mat(c_0+i,j)=ic
*		  ic=c_1
*	       endif
	       ic=ic+1
	    enddo
	 enddo
      endif
600   format ('... Quantisizing to ',i3,' x ',i3,' shades ...')
      end
