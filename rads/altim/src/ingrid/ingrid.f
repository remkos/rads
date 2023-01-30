      program ingrid
***********************************************************************
* The ultimate grid facility :                                        *
*                                                                     *
*                         I n G r i d                                 *
*                                                                     *
* Remko Scharroo, Delft, 15 September 1994                            *
***********************************************************************
      include "ingrid.inc"
      integer mfunc
      parameter (mfunc = 8+2*MPOL)
      real*4 funpar(mfunc),fi(MBUF)
      real*8 zmean,zrms,sumx,sumxx,sumy,sumyy,sumxy
      logical error,end5,quit,stdfmt,fresh6/.false./,sg,change,listen
      integer ifunpa(4),kbi(MBUF),idi(MBUF)
      character*80 answer,filenm,fnmenu,fntemp
      character optio(MOPT,MOPT)*41
      data funpar/1.,1.,1.,1.,0.,1.,0.,1.,
     .1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0./,
     .ifunpa/0,0,0,0/
      parameter (fntemp='ingrid.tmp')

      integer kx,ky,nxp,nyp,nxc,nyc,nxt,nyt,mx,my,lox,loy,
     .minx,miny,newx,newy,maxx,maxy,nx1,ny1,nx2,ny2
      integer itype,igeo,mxy,ipos,nroom,itest,imax,iopt,istep,
     .ndeg,ndeg1,ndeg2,nspace,ni,indexr,indexi,number,
     .menu,modifd,lnblnk
      integer i,j,k,m,l
      integer i2,i3,id1,id2,id3/0/,id4,kb,kb1,kb2,kb3,kb4
      real f1,f2,f3,b0,bx,by,bxy,fx,fy,ff1,ff2,frx,fry,bili
      real x,y,z,expand,pi,val,angle,atemp,fact,tilt,bias,
     .wieber,x0,x1,y0,y1
      real zmin,zmax,corr(MBUF,MBUF)
*
  550 format (a)
  551 format (a,$)
*
* Initialization
*
      pi=4*atan(1.)
      do 10 kb=1,MBUF
   10    used(kb)=.false.
      menu=1
      itest=1
      stdfmt=.true.
      listen=.true.
      end5=.true.
*
* Read menu headers and options.
* optio(j,i): i=menu#; j=1:menuhead; j=option#+2
*
      fnmenu='/user/altim'
      call checkenv('ALTIM',fnmenu,l)
      open(10,file=fnmenu(:l)//'/src/ingrid/ingrid.menus',
     |status='old',err=1350)
      rewind(10)
      do 11 i=1,MOPT
         do 11 j=1,MOPT
   11       read(10,550,err=12,end=12) optio(j,i)
*
   12 idnext=1
      close (10)
      write (0,1010) MVRS
*
* List buffer information
*
 1010 format ('This is InGrid version ',i4.4)
 1019 format (/'Buf |  x-range  | y-range | dimension | file',$)
 1020 format (/i3,' |',2i5,' |',2i4,' |',2i5,' |',$)
 1030 format (/'None of the ',i3,' buffers is used',$)
 1040 format (/'Number of buffers used :',i3,$)
 1050 format (' *NEW*',$)
 1060 format (' ',a,$)
 1070 format (' (used)',$)
*
   20 if (itest.eq.0) menu=1
      itest=1
      jused=0
      idnext=1
      number=0
      do 25 kb=1,MBUF
         blkd(kb)=.false.
         if (used(kb)) then
            number=number+1
            jused=jused+nx(kb)*ny(kb)
            idnext=max(idnext,id(kb)+nx(kb)*ny(kb))
            if (listen) then
               if (number.eq.1) write (0,1019)
               write (0,1020) kb,nint(xmin(kb)),nint(xmax(kb)),
     1		nint(ymin(kb)),nint(ymax(kb)),nx(kb),ny(kb)
               if (fname(kb).eq.' ') then
                  write (0,1050)
               else
		  filenm=fname(kb)
		  l=lnblnk(filenm)
                  write (0,1060) filenm(:l)
               endif
               if (proc(kb)) write (0,1070)
            endif
         endif
   25 continue
      if (number.eq.0) then
         write (0,1030) MBUF
      else if (.not.listen) then
         write (0,1040) number
      endif
      if (jused.ne.idnext-1) then
	 write (0,550)
	 write (0,550)
	 write (0,550) 'Repacking buffers ...'
	 call repack
      endif
      write (0,550)
      goto (100,200,300,400,500,530,600,700),menu
*
* Print options menu #1, which is main menu
*
  100 call menuh(iopt,optio(1,1),quit)
      if (quit) goto 190
      menu=iopt
      fresh6=.true.
*
* Options 1 and 2 go to same menu (#2)
*
      if (iopt.eq.1) menu=2
      if (menu.eq.2) goto 200
      if (menu.ge.3) goto 20
*
* Here clean exit (after choosing option 0 in main menu).
* Label 190 is called at all 'end=' switches for user input.
*
 1192 format ('ingrid:',i3,' unwritten buffer(s) lossed')
 1194 format ('ingrid:',i3,' unwritten buffer(s)'/
     .8x,'OK to lose changes? -> ',$)
*
      end5=.false.
  190 continue
      if (end5) write (0,550) 'ingrid: EOD on standard input'
      modifd=0
      do 192 kb=1,MBUF
  192    if (used(kb) .and. fname(kb).eq.' ' .and.
     1		.not.proc(kb)) modifd=modifd+1
      if (modifd.gt.0) then
         if (end5) then
            write (0,1192) modifd
         else
            end5=.true.
            write (0,1194) modifd
            read (5,550,end=190,err=194) answer
            if (answer(1:1).eq.'n') goto 20
         endif
      endif
  194 close (11,status='delete')
      goto 999
*
* Option 2: Read & write submenu
*
  200 continue
      menu=1
      if (stdfmt) then
         itest=1
      else
         call menuh(itest,optio(1,2),quit)
         if (quit) goto 190
         if (itest.eq.0) goto 20
      endif
      if (iopt.eq.2) goto 240
*
* Here reading file(s) in one of three formats
*
      write (0,551) 'Enter file name (wildcards allowed) -> '
      read (5,550,end=190,err=20) filenm
*
* Make list of candidate files
*
      call system('ls '//filenm//'>'//fntemp)
      open (11,file=fntemp,status='old',err=230)
      rewind (11)
*
* Open candidate files and read
*
      call ireadf(itest,quit)
      if (quit) goto 190
  230 close (11,status='delete')
      goto 20
*
* Here writing file in one of three formats
*
  240 continue
      call oldbuf('Enter buffer nr.',kb,.false.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      write (0,551) 'Enter file name -> '
      read (5,550,end=190,err=20) filenm
      call iwritef(filenm,itest,kb,nx(kb),ny(kb),quit)
      if (quit) goto 190
      close (10)
      goto 20
*
* Option 3: Buffer manipulations submenu
*
  300 continue
      call menuh(itest,optio(1,3),quit)
      if (quit) goto 190
      goto (20,310,320,330,340,350,360,370),itest+1
*
* Copy buffer
*
  310 continue
      call twobuf(kb1,kb2,.false.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      id1=id(kb1)-1
      id2=id(kb2)-1
      nspace=nx(kb1)*ny(kb1)
      do 315 i=1,nspace
  315    a(id2+i)=a(id1+i)
      proc(kb1)=.true.
      fname(kb2)=' '
      proc(kb2)=.false.
      used(kb2)=.true.
      call cparea(kb1,kb2)
      goto 20
*
* Move buffer
*
  320 continue
      call oldbuf('Enter source buffer',kb1,.false.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      write (0,551) 'Enter target buffer -> '
      read (5,*,err=20,end=190) kb2
      if (kb1.eq.kb2) then
         write (0,550)
     .   'ingrid: No use for buffer to be both source and target'
         goto 20
      endif
      call noclob(kb2,error,quit)
      if (quit) goto 190
      if (error) goto 20
      id(kb2)=id(kb1)
      nx(kb2)=nx(kb1)
      ny(kb2)=ny(kb1)
      used(kb2)=used(kb1)
      fname(kb2)=' '
      proc(kb2)=proc(kb1)
      fname(kb2)=fname(kb1)
      call cparea(kb1,kb2)
      used(kb1)=.false.
      goto 20
*
* Clear buffer
*
  330 continue
 1330 format ('Buffer',i3,': Cleared')
      call multin('Enter buffer(s) to be cleared',MBUF,ni,kbi,error,
     .quit)
      if (quit) goto 190
      if (error) goto 20
      do 335 i=1,ni
         kb=iabs(kbi(i))
         call noclob(kb,error,quit)
         if (quit) goto 190
         if (kb.eq.0) goto 20
         if (error) goto 335
         if (used(kb)) write (0,1330) kb
         used(kb)=.false.
  335 continue
      goto 20
*
* Mark all buffers as unchanged
*
  340 continue
      do kb=1,MBUF
         if (used(kb) .and. fname(kb).eq.' ')
     1      fname(kb)='(unchanged)'
      enddo
      goto 20
*
* Save changed files to /tmp
*
  350 continue
      write (0,550) 'ingrid: option removed'
      goto 20
*
* Re-pack buffers
*
  360 continue
      call repack
      goto 20
*
* Show buffer status
*
  370 continue
      number=0
      do 371 kb=1,MBUF
         if (used(kb)) then
            number=number+1
            kbi(number)=kb
         endif
  371 continue
      call bubble(kbi,id,number)
      if (number.gt.0) write (0,1370)
      do 375 k=1,number
         kb=kbi(k)
         nspace=nx(kb)*ny(kb)
         if (k.eq.number) then
            nroom=idnext-id(kb)
         else
            nroom=id(kbi(k+1))-id(kb)
         endif
  375    write (0,1372) kb,nx(kb),ny(kb),nspace,nroom,
     .int(100./MA*nspace+.5),int(100./MA*nroom+.5),
     .int(100./nroom*nspace+.5),id(kb)
      nroom=idnext-1
      if (number.gt.0) then
         write (0,1374) jused,nroom,int(100./MA*jused+.5),
     .   int(100./MA*nroom+.5),int(100./nroom*jused+.5),MA
      endif
      if (number.gt.5) then
         write (0,"(/'... Press Enter-key to continue ...')")
         read (5,550,err=20,end=190) answer
      endif
 1370 format (/'Buffer',t15,'grid',t30,'nr. of points',t47,'% of total',
     .t59,'occupation',t73,'start'/'  nr.',t13,'dimension',t32,
     .'grid-buffer',t47,'grid-buff',t62,'rate',t72,'address')
 1372 format (1x,i3,t10,'(',i6,',',i6,')',t29,i7,'-',i7,
     .t48,i3,'-',i3,'%',t62,i3,'%',t72,i7)
 1374 format ('total',t29,i7,'-',i7,
     .t48,i3,'-',i3,'%',t62,i3,'%',t72,i7)
      goto 20
*
* Option 4: Grid size manipulations submenu
*
  400 continue
      call menuh(itest,optio(1,4),quit)
      if (quit) goto 190
      goto (20,410,420,420,430,440,460),itest+1
*
* Select subsection of grid
*
  410 continue
      kb2=0
      call oldbuf('Enter source buffer',kb1,.true.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      if (kb1.lt.0) then
         kb1=-kb1
         kb2=kb1
      endif

      write (0,551) 'Enter x-range or x-index-range of subsection -> '
      nx1=nx(kb1)
      call getrange(minx,maxx,x0,x1,
     1	xmin(kb1),xmax(kb1),nx1,quit,error)
      if (quit) goto 190
      if (error) goto 20
      write (0,416) 'x',minx,maxx,x0,x1

      write (0,551) 'Enter y-range or y-index-range of subsection -> '
      ny1=ny(kb1)
      call getrange(miny,maxy,y0,y1,
     1	ymin(kb1),ymax(kb1),ny1,quit,error)
      if (quit) goto 190
      if (error) goto 20
      write (0,416) 'y',miny,maxy,y0,y1

      lox=minx-2
      loy=miny-2
      newx=1+maxx-minx
      newy=1+maxy-miny
      if (kb2.eq.0) then
         call newbuf('Enter target buffer',kb2,newx,newy,
     .   error,quit)
         if (quit) goto 190
         if (error) goto 20
      endif
      do 415 ky=1,newy
        do 415 kx=1,newx
  415     a(id(kb2)+(kx-1)+(ky-1)*newx)=
     .       a(id(kb1)+(kx+lox)+(ky+loy)*nx1)
      proc(kb1)=.true.
      fname(kb2)=' '
      proc(kb2)=.false.
      used(kb2)=.true.
      xmin(kb2)=x0
      xmax(kb2)=x1
      ymin(kb2)=y0
      ymax(kb2)=y1
416   format (a1,': index :',i4,i5,'  range : ',2f8.2)
      goto 20
*
* Expand grid by integer number or change grid dimensions
* Combination of sub-options 2 and 3
*
  420 continue
      call oldbuf('Enter source buffer',kb1,.false.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      nx1=nx(kb1)
      ny1=ny(kb1)
      if (itest.eq.2) then
         write (0,551) 'Enter expansion factor -> '
         read (5,*,end=190,err=20) expand
         if (expand.lt.0) then
            write (0,550) 'ingrid: illegal expansion factor'
            goto 20
         endif
         nx2=int(expand*(nx1-1)+1.5)
         ny2=int(expand*(ny1-1)+1.5)
      else
         write (0,551) 'Enter new x- and y-dimensions -> '
         read (5,*,end=190,err=20) nx2,ny2
         if (nx2.lt.1 .or. ny2.lt.1) then
            write (0,550) 'ingrid: at least one dimension is illegal'
            goto 20
         endif
      endif
      blkd(kb1)=.true.
      call newbuf('Enter target buffer',kb2,nx2,ny2,error,quit)
      if (quit) goto 190
      if (error) goto 20
      call grsexp(a(id(kb1)),a(id(kb2)),nx1,ny1,nx2,ny2,nx1,nx2,error)
      if (error) then
         write (0,550) 'ingrid: target buffer lost'
         used(kb2)=.false.
      else
         proc(kb1)=.true.
         fname(kb2)=' '
         proc(kb2)=.false.
         used(kb2)=.true.
      endif
      call cparea(kb1,kb2)
      goto 20
*
* Make a patchwork of a number of grids
*
  430 continue
 1430 format ('ingrid: too few (< 0) or too many (>',i3,') patches')
      change=.false.
      nxp=0
      nyp=0
      write (0,551) 'Enter number of grids in horizontal direction -> '
      read (5,*,end=190,err=20) mx
      write (0,551) 'Enter number of grids in vertical direction   -> '
      read (5,*,end=190,err=20) my
      if (mx.le.0 .or. my.le.0 .or. mx*my.gt.MBUF) then
         write (0,1430) MBUF
         goto 20
      endif
      call multin('Enter buffer numbers, left to right, top to bottom'//
     .' (0 for empty)',MBUF,mxy,kbi,error,quit)
      if (quit) goto 190
      if (error) goto 20
      write (0,"(3x,'configuration: ',$)")
      do 432 ky=my,1,-1
         if (ky.lt.my) write (0,"(/18(' '),$)")
         do 432 kx=1,mx
            i=kx+(my-ky)*mx
            if (.not.used(kbi(i)) .or. i.gt.mxy) kbi(i)=0
            kb=kbi(i)
            if (kb.eq.0) then
               change=.true.
            else
               blkd(kb)=.true.
               if (nxp*nyp.eq.0) then
                  nxp=nx(kb)
                  nyp=ny(kb)
               else if (nx(kb).ne.nxp .or. ny(kb).ne.nyp) then
                  write (0,550)
     .            'ingrid: dimensions must be equal for all grids'
                  goto 20
               endif
            endif
  432       write (0,"(i3,$)") kb
      write (0,550)
      if (nxp*nyp.eq.0) then
         write (0,550) 'ingrid: no used grids, action terminated'
         goto 20
      endif
      call newbuf('Enter buffer number for result',kb1,
     .nxp*mx,nyp*my,error,quit)
      if (quit) goto 190
      if (error) goto 20
      id1=id(kb1)-1
      if (change) then
         write (0,551)
     .   'Enter value to be substituted in empty patches -> '
         read (5,*,end=190,err=20) x
         do 436 i=1,nxp*nyp*mx*my
  436       a(id1+i)=x
      endif
      do 438 kx=1,mx
         do 438 ky=my,1,-1
            kb=kbi(kx+(my-ky)*mx)
            if (kb.ne.0) then
               call subarr(a(id(kb)),nxp,nyp,
     .         a(id(kb1)),nx(kb1),(kx-1)*nxp+1,(ky-1)*nyp+1)
               proc(kb)=.true.
            endif
  438 continue
      fname(kb1)=' '
      proc(kb1)=.false.
      used(kb1)=.true.
      goto 20
*
* Overlapping
*
  440 continue
 1440 format (/
     .'  A°°°°°°°       A°°°±°°°           °°°°       A°°°°    '/
     .'   ±±±±±±±        °°°±°°°       A°°°±°°°B       °°°±°°° '/
     .'   °°°°°°°B       °°°±°°°B       °°°°              °°°°B'/
     .'      1              2              3              4    '/
     .'Select proper overlap geometry (1-4) -> ',$)
      write (0,1440)
      read (5,*,end=190,err=20) igeo
      if (igeo.lt.1 .or. igeo.gt.4) goto 20
      call oldbuf('Enter source buffer A',kb1,.false.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      nx1=nx(kb1)
      ny1=ny(kb1)
      call oldbuf('Enter source buffer B',kb2,.false.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      nx2=nx(kb2)
      ny2=ny(kb2)
      if (igeo.eq.1) then
         if (nx1.ne.nx2) then
            write (0,550)
     .'ingrid: horizontal (x)dimension of grids A and B must be equal'
            goto 20
         endif
         nxt=nx1
         nxc=nx1
      else
         write (0,551) 'Enter size of overlap in x-direction'//
     .   ' (horiz.) (number of cells) -> '
         read (5,*,end=190,err=20) nxc
         if (nxc.gt.nx1 .or. nxc.gt.nx2) then
            write (0,550) 'ingrid: horizontal overlap too large'
            goto 20
         endif
         nxt=nx1+nx2-1-nxc
      endif
      if (igeo.eq.2) then
         if (ny1.ne.ny2) then
            write (0,550)
     .'ingrid: vertical (y)dimension of grids A and B must be equal'
            goto 20
         endif
         nyt=ny1
         nyc=ny1
      else
         write (0,551) 'Enter size of overlap in y-direction'//
     .   ' (verti.) (number of cells) -> '
         read (5,*,end=190,err=20) nyc
         if (nyc.gt.ny1 .or. nyc.gt.ny2) then
            write (0,550) 'ingrid: vertical overlap too large'
            goto 20
         endif
         nyt=ny1+ny2-1-nyc
      endif
      if (igeo.ge.3) then
         write (0,551)
     .   'Which value must be substituted in non-covered areas -> '
         read (5,*,end=190,err=20) x
      endif
      itype=1
      if (nxc.ge.0 .and. nyc.ge.0) then
         write (0,550) 'Enter type of overlap'
         write (0,551)
     .   '1. A over B   2. B over A   3. average   4. smooth   -> '
         read (5,*,end=190,err=20) itype
      endif
      if ((nxc.eq.0 .or. nyc.eq.0) .and. itype.eq.4) itype=3
      if (itype.lt.1 .or. itype.gt.4) then
         write (0,550) 'ingrid: illegal entry'
         goto 20
      endif
      blkd(kb1)=.true.
      blkd(kb2)=.true.
      call newbuf('Enter target buffer',kb3,nxt,nyt,error,quit)
      if (quit) goto 190
      if (error) goto 20
*
      nspace=nxt*nyt-1
      id1=id(kb1)
      id2=id(kb2)
      id3=id(kb3)
*
* Substitute background value
*
*     write (6,*) igeo,id3,nspace,x
      if (igeo.ge.3) then
*	 write (6,*) 'Filling space with',x
         do 442 i=0,nspace
  442       a(id3+i)=x
      endif
      if (itype.eq.2) then
*
* Substitute first A then B
*
         my=1
         if (igeo.eq.1 .or. igeo.eq.4) my=nyt-(ny1-1)
         call subarr(a(id1),nx1,ny1,a(id3),nxt,1,my)
         mx=1
         my=1
         if (igeo.ge.2) mx=nxt-(nx2-1)
         if (igeo.eq.3) my=nyt-(ny2-1)
         call subarr(a(id2),nx2,ny2,a(id3),nxt,mx,my)
      else
*
* Substitute first B then A
*
         mx=1
         my=1
         if (igeo.ge.2) mx=nxt-(nx2-1)
         if (igeo.eq.3) my=nyt-(ny2-1)
         call subarr(a(id2),nx2,ny2,a(id3),nxt,mx,my)
         my=1
         if (igeo.eq.1 .or. igeo.eq.4) my=nyt-(ny1-1)
         call subarr(a(id1),nx1,ny1,a(id3),nxt,1,my)
      endif
      if (itype.ge.3) then
*
* Make smooth overlap or take the average
*
         if (igeo.eq.1) then
            do 451 ky=0,nyc
               if (itype.eq.3) then
                  f2=.5
               else
                  f2=(1+cos(ky*pi/nyc))/2
               endif
               do 451 kx=1,nxt
                  i=(nyt-ny1+ky)*nxt+kx-1
  451             a(id3+i)=a(id2+i)*f2+a(id3+i)*(1-f2)
         else if (igeo.eq.2) then
            do 452 kx=0,nxc
               if (itype.eq.3) then
                  f3=.5
               else
                  f3=(1+cos(kx*pi/nxc))/2
               endif
               do 452 ky=1,nyt
                  i2=(ky-1)*nx2+kx+id2
                  i3=(ky-1)*nxt+(nxt-nx2+kx)+id3
  452             a(i3)=a(i3)*f3+a(i2)*(1-f3)
         else if (igeo.eq.3) then
            do 453 kx=0,nxc
               do 453 ky=0,nyc
                  if (itype.eq.3) then
                     f3=.5
                  else
                     f3=wieber(kx,nxc,ky,nyc,pi)
                  endif
                  i2=kx+ky*nx2+id2
                  i3=(nyt-ny2+ky)*nxt+(nxt-nx2+kx)+id3
  453             a(i3)=a(i3)*f3+a(i2)*(1-f3)
         else
            do 454 kx=0,nxc
               do 454 ky=0,nyc
                  if (itype.eq.3) then
                     f3=.5
                  else
                     f3=wieber(kx,nxc,nyc-ky,nyc,pi)
                  endif
                  i2=kx+(nyt-ny1+ky)*nx2+id2
                  i3=(nyt-ny1+ky)*nxt+(nxt-nx2+kx)+id3
  454             a(i3)=a(i3)*f3+a(i2)*(1-f3)
         endif
      endif
      proc(kb1)=.true.
      proc(kb2)=.true.
      fname(kb3)=' '
      proc(kb3)=.false.
      used(kb3)=.true.
      goto 20
*
* Rotate grid
*
  460 continue
      call twobuf(kb1,kb2,.false.,error,quit)
      write (0,"('Rotate grid...   1.  90 deg clockwise'/
     .17x,'2.  90 deg counter-clockwise'/
     .17x,'3. 180 degrees               -> ',$)")
      read (5,*,end=190,err=20) ipos
      if (ipos.lt.1 .or. ipos.gt.3) then
         write (0,550) 'ingrid: illegal entry'
         goto 20
      else if (ipos.eq.3) then
         nspace=nx(kb1)*ny(kb1)
         id1=id(kb1)-1
         id2=id(kb2)+nspace
         do 462 i=1,nspace
  462       a(id2-i)=a(id1+i)
      else
         call rot90(ipos,a(id(kb1)),a(id(kb2)),nx(kb1),ny(kb1))
         nx(kb2)=ny(kb1)
         ny(kb2)=nx(kb1)
      endif
      proc(kb1)=.true.
      fname(kb2)=' '
      proc(kb2)=.false.
      used(kb2)=.true.
      goto 20
*
* Option 5: Data manipulations on one grid
*
  500 continue
      call menuh(itest,optio(1,5),quit)
      if (quit) goto 190
      goto (20,710,510,510,510,510,510,510,510,510),itest+1
*
* Transformation of one grid
*
  510 continue
*
* Combination of sub-options 2, 3, 4, 5, 6, 7, and 8.
* Transformations of one grid.
* 2 = linear transformation of grid (tilt/bias)
* 3 = linear transformation of grid (tilt/bias) using mask
* 4 = compute square of grid
* 5 = take square root of grid
* 6 = subtract mean from grid
* 7 = illuminate grid
* 8 = reciprocal of grid
* 9 = 10log of grid
*
      call twobuf(kb1,kb2,(itest.ne.7),error,quit)
      if (quit) goto 190
      if (error) goto 20
      zmean=0
      nxp=nx(kb1)
      nyp=ny(kb1)
      nspace=nxp*nyp
      id1=id(kb1)-1
      id2=id(kb2)-1
      if (itest.eq.3) then
         write (0,550) 'Applied only to values between two limits'
         write (0,551) 'Enter low and high limits -> '
         read (5,*,end=190,err=20) zmin,zmax
         if (zmin.gt.zmax) then
            write (0,550) 'ingrid: high limit less than low limit'
            goto 20
         endif
      endif
      if (itest.le.3) then
         write (0,551) 'Enter offset -> '
         read (5,*,end=190,err=20) bias
         write (0,551) 'Enter factor -> '
         read (5,*,end=190,err=20) tilt
      endif
      if (itest.eq.2) then
         do 520 i=1,nspace
  520       a(id2+i)=bias+tilt*a(id1+i)
      else if (itest.eq.3) then
         do 522 i=1,nspace
            val=a(id1+i)
            if (val.ge.zmin .and. val.le.zmax) val=bias+tilt*val
  522       a(id2+i)=val
      else if (itest.eq.4) then
         do 524 i=1,nspace
  524       a(id2+i)=a(id1+i)**2
      else if (itest.eq.5) then
         do 526 i=1,nspace
	    val=a(id1+i)
	    if (val.lt.0 .or. val.gt.1e20) then
	       a(id2+i)=1e35
            else
	       a(id2+i)=sqrt(val)
	    endif
  526    continue
      else if (itest.eq.6) then
         do 527 i=1,nspace
  527       zmean=zmean+a(id1+i)
         zmean=zmean/nspace
         do 528 i=1,nspace
  528       a(id2+i)=a(id1+i)-zmean
      else if (itest.eq.7) then
         if (nxp.lt.3 .or. nyp.lt.3) then
            write (0,550) 'ingrid: grid to small to be illuminated'
            goto 20
         endif
         write (0,551) 'Enter illumination direction'//
     .   ' (deg, <- = 0, counter-clockwise) -> '
         read (5,*,err=20,end=190) angle
         call grillu(a(id1+1),a(id2+1),angle,nxp,nyp,nxp)
      else if (itest.eq.8) then
	 do 529 i=1,nspace
	    val=a(id1+i)
            if (val.gt.1e20 .or. val.lt.1e-20) then
	       a(id2+i)=1e35
	    else
	       a(id2+i)=1/val
	    endif
  529    continue
      else if (itest.eq.9) then
         do i=1,nspace
	    val=a(id1+i)
	    if (val.le.0 .or. val.gt.1e20) then
	       a(id2+i)=1e35
	    else
	       a(id2+i)=log10(val)
	    endif
	 enddo
      endif
      proc(kb1)=.true.
      fname(kb2)=' '
      proc(kb2)=.false.
      used(kb2)=.true.
      call cparea(kb1,kb2)
      goto 20
*
* Option 6: Linear combination of several grids
*
  530 continue
      call menuh(itest,optio(1,6),quit)
      if (quit) goto 190
      goto (20,710,540,540,540,560,570,580,540),itest+1
*
* Transformations with several grids
*
  540 continue
*
* Combination of sub-options 2, 3, 4, 8, and 9
* 2 = linear combination of several grids
* 3 = mean grid of several grids
* 4 = rms grid of several grids
* 8 = correlation between several grids
* 9 = product of several grids
*
      call multin('Enter source buffer(s)',MBUF,ni,kbi,error,quit)
      if (quit) goto 190
      if (error) goto 20
      nxt=nx(kbi(1))
      nyt=ny(kbi(1))
      nspace=nxt*nyt
      do 542 i=1,ni
         kb=kbi(i)
         if (kb.lt.1 .or. kb.gt.MBUF) then
            write (0,550) 'ingrid: at least one illegal buffer number'
            goto 20
         else if (.not.used(kb)) then
            write (0,550) 'ingrid: at least one buffer not in use'
            goto 20
         else if (nx(kb).ne.nxt .or. ny(kb).ne.nyt) then
            write (0,550)
     .      'ingrid: error; buffers must have equal dimensions'
            goto 20
         endif
         fi(i)=1./ni
         idi(i)=id(kb)-1
  542 continue
      if (itest.eq.2) then
         write(0,"('Enter ',i2,' factors -> ',$)") ni
         read(5,*,end=190,err=20) (fi(i),i=1,ni)
      endif
      if (itest.eq.8) then
         call newbuf('Enter target buffer for correlation grid',
     .   kb3,ni,ni,error,quit)
         if (.not.error) id3=id(kb3)-1
      else
         call newbuf('Enter target buffer for result',
     .   kb3,nxt,nyt,error,quit)
         if (error) goto 20
         id3=id(kb3)-1
      endif
      if (quit) goto 190
      if (itest.eq.4) then
         do 546 j=1,nspace
            atemp=0
            do 544 i=1,ni
  544          atemp=atemp+a(idi(i)+j)**2
  546       a(id3+j)=sqrt(atemp/ni)
      else if (itest.eq.8) then
         do i=1,ni
            do k=1,i
	       sumx=0d0
	       sumy=0d0
	       sumxx=0d0
	       sumxy=0d0
	       sumyy=0d0
	       l=0
               do j=1,nspace
		  x=a(idi(i)+j)
		  y=a(idi(k)+j)
		  if (abs(x)+abs(y).lt.1d20) then
		     sumx=sumx+x
		     sumy=sumy+y
		     sumxx=sumxx+x*x
		     sumxy=sumxy+x*y
		     sumyy=sumyy+y*y
		     l=l+1
		  endif
	       enddo
	       if (l.gt.0) then
	          sumx=sumx/l
	          sumy=sumy/l
	          corr(i,k)=(sumxy/l-sumx*sumy)/
     &			sqrt(sumxx/l-sumx*sumx)/
     &			sqrt(sumyy/l-sumy*sumy)
	       else
		  corr(i,k)=1e30
	       endif
	       corr(k,i)=corr(i,k)
	    enddo
	 enddo
         istep=18
         write (6,599)
         do 596 i=1,ni,istep
            imax=min0(i+istep,ni)-1
            do 594 k=i+1,ni
               write (6,"(/i3,'|',$)") kbi(k)
               do 594 m=i,min0(imax,k-1)
                  if (.not.error) then
                     a(id3+k+(m-1)*ni)=corr(m,k)
                     a(id3+m+(k-1)*ni)=0
                  endif
  594             write (6,"(i4,$)") nint(corr(m,k)*100)
            write (6,"(/'   +',18(a,:))") ('----',j=i,imax)
  596       write (6,"(4x,18(i4,:))") (kbi(j),j=i,imax)
      else if (itest.eq.9) then
         do 598 j=1,nspace
            atemp=1
            do 597 i=1,ni
	       val=a(idi(i)+j)
	       if (val.gt.1e20) then
		  atemp=1e35
		  goto 598
	       endif
  597          atemp=atemp*val
  598       a(id3+j)=atemp
      else
         do 554 j=1,nspace
            atemp=0
            do 552 i=1,ni
	       val=a(idi(i)+j)
	       if (val.gt.1e20) then
		  atemp=1e35
		  goto 554
	       endif
  552          atemp=atemp+fi(i)*val
  554       a(id3+j)=atemp
      endif
      do 556 i=1,ni
  556    proc(kbi(i))=.true.
      if (.not.error) then
         fname(kb3)=' '
         proc(kb3)=.false.
         used(kb3)=.true.
	 call cparea(kbi(1),kb3)
      endif
  599 format(/'Correlations (%)')
      goto 20
*
* Difference of several grids from reference grid
*
  560 continue
      call oldbuf('Enter buffer number of reference grid',
     .kb3,.false.,error,quit)
      blkd(kb3)=.true.
      id3=id(kb3)-1
      write (0,550) 'Note: the source buffer(s) that you will be',
     .'entering beneath are also target buffer(s)'
      call multin('Enter source/target buffer(s)',MBUF,ni,kbi,error,
     .quit)
      if (quit) goto 190
      if (error) goto 20
      nxt=nx(kbi(1))
      nyt=ny(kbi(1))
      nspace=nxt*nyt
      change=.false.
      do 562 i=1,ni
         kb=kbi(i)
         if (kb.lt.1 .or. kb.gt.MBUF) then
            write (0,550) 'ingrid: at least one illegal buffer number'
            goto 20
         else if (blkd(kb)) then
            write (0,"('ingrid: you cannot reuse buffer ',i3)") kb
            goto 20
         else if (.not.used(kb)) then
            write (0,550) 'ingrid: at least one buffer not in use'
            goto 20
         else if (nx(kb).ne.nxt .or. ny(kb).ne.nyt) then
            write (0,550)
     .      'ingrid: error; buffers must have equal dimensions'
            goto 20
         endif
         if (fname(kb).eq.' ' .and. .not.proc(kb)) change=.true.
         idi(i)=id(kb)-1
  562 continue
      if (change) then
         write (0,550) 'ingrid: Unwritten changes will be lost'
         write (0,551) '        Continue? -> '
         read (5,550,end=190,err=20) answer
         if (answer(1:1).eq.'n') goto 20
      endif
      do 564 i=1,ni
         fname(kbi(i))=' '
         proc(kbi(i))=.false.
         do 564 j=1,nspace
  564       a(idi(i)+j)=a(idi(i)+j)-a(id3+j)
      proc(kb3)=.true.
      call cparea(kbi(1),kb3)
      goto 20
*
* Fit one grid with another
*
  570 continue
      call oldbuf('Enter buffer to be approximated',
     .kb1,.true.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      kb3=0
      if (kb1.lt.0) then
         kb1=-kb1
         kb3=kb1
      endif
      call oldbuf('Enter buffer containing approximation function',
     .kb2,.false.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      nxt=nx(kb1)
      nyt=ny(kb1)
      nspace=nxt*nyt
      if (nxt.ne.nx(kb2) .or. nyt.ne.ny(kb2)) then
         write (0,550)
     .   'ingrid: error; grids must have equal dimensions'
         goto 20
      endif
      if (kb3.eq.0) then
         call newbuf('Enter target buffer for result',
     .   kb3,nxt,nyt,error,quit)
         if (quit) goto 190
         if (error) goto 20
      endif
      call grest1(a(id(kb1)),a(id(kb2)),fact,nxt,nyt,nxt)
      write (0,1571) fact
 1571 format ('Estimated coefficient: ',f12.4)
      id2=id(kb2)-1
      id3=id(kb3)-1
      do 575 i=1,nspace
  575    a(id3+i)=a(id2+i)*fact
      proc(kb1)=.true.
      proc(kb2)=.true.
      fname(kb3)=' '
      proc(kb3)=.false.
      used(kb3)=.true.
      call cparea(kb1,kb3)
      goto 20
*
* Estimate best fitting linear combination of two other grids
*
  580 continue
      call oldbuf('Enter buffer to be approximated',
     .kb1,.true.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      kb4=0
      if (kb1.lt.0) then
         kb1=-kb1
         kb4=kb1
      endif
      write (0,551)
     .'Enter two buffers containing approximation functions -> '
      read (5,*,end=190,err=20) kb2,kb3
      if (min(kb2,kb3).lt.1 .or. max(kb2,kb3).gt.MBUF) then
         write (0,550) 'ingrid: at least one illegal buffer number'
         goto 20
      endif
      if (.not.(used(kb2) .and. used(kb3))) then
         write (0,550) 'ingrid: at least one buffer not in use'
         goto 20
      endif
      nxp=nx(kb1)
      nyp=ny(kb1)
      nspace=nxp*nyp
      if (nxp.ne.nx(kb2) .or. nxp.ne.nx(kb3) .or.
     .   nyp.ne.ny(kb2) .or. nyp.ne.ny(kb3)) then
         write (0,550)
     .   'ingrid: error; grids must have equal dimensions'
         goto 20
      endif
      call grest2(a(id(kb1)),a(id(kb2)),f2,
     .a(id(kb3)),f3,nxp,nyp,nxp,sg)
      if (sg) then
         write (0,550) 'ingrid: singular system of equations'
         goto 20
      endif
      write (0,1581) f2,f3
 1581 format ('Estimated respective coefficients: ',2f12.4)
      if (kb4.eq.0) then
         call newbuf('Enter target buffer for result',kb4,nxp,nyp,
     .   error,quit)
         if (quit) goto 190
         if (error) goto 20
      endif
      id2=id(kb2)-1
      id3=id(kb3)-1
      id4=id(kb4)-1
      do 585 i=1,nspace
  585    a(id4+i)=a(id2+i)*f2+a(id3+i)*f3
      proc(kb1)=.true.
      proc(kb2)=.true.
      proc(kb3)=.true.
      fname(kb4)=' '
      proc(kb4)=.false.
      used(kb4)=.true.
      call cparea(kb1,kb4)
      goto 20
*
* Option 7: Create grid of function
*
  600 continue
      if (fresh6) write (0,1603)
      fresh6=.false.
 1603 format (/'------Function will be-----------'/
     .'   B(kx,ky) * F1(f1(kx)) * F2(f2(ky))     where:'/
     .' - kx and ky are the grid point numbers'/
     .' - f1 and f2 are linear relations, both default 0 to 1'/
     .' - B is bilinear,                       default 1'/
     .' - F1 and F2 are other functions,  both default 1')
      call menuh(itest,optio(1,7),quit)
      if (quit) goto 190
      goto (20,610,620,630,640,650,660,670,680,690),itest+1
*
* Change bilinear
*
  610 continue
      write (0,550) 'Enter function values at corner points of grid'
      write (0,551) 'Left and right for low  y (2 values) -> '
      read (5,*,end=190,err=20) funpar(1),funpar(2)
      write (0,551) 'Left and right for high y (2 values) -> '
      read (5,*,end=190,err=20) funpar(3),funpar(4)
      goto 20
*
* Change f1 low and high values
*
  620 continue
      write (0,551) 'Enter min. and max. of f1 (2 values) -> '
      read (5,*,end=190,err=20) funpar(5),funpar(6)
      goto 20
*
* Change f2 low and high values
*
  630 continue
      write (0,551) 'Enter min. and max. of f2 (2 values) -> '
      read (5,*,end=190,err=20) funpar(7),funpar(8)
      goto 20
*
* Toggle sin(pi * f1) application
*
  640 continue
      if (ifunpa(1).eq.0) then
         write (0,550) 'sin(pi * f1) will be used instead of f1'
         ifunpa(1)=1
      else
         write (0,550) 'sine computation for f1 disabled'
         ifunpa(1)=0
      endif
      goto 20
*
* Toggle sin(pi * f2) application
*
  650 continue
      if (ifunpa(2).eq.0) then
         write (0,550) 'sin(pi * f2) will be used instead of f2'
         ifunpa(2)=1
      else
         write (0,550) 'sine computation for f2 disabled'
         ifunpa(2)=0
      endif
      goto 20
*
* Change F1 polynomial
*
  660 continue
      indexi=3
      indexr=9
      goto 671
*
* Change F2 polynomial
*
  670 continue
      indexi=4
      indexr=16
  671 ifunpa(indexi)=0
      write (0,551) 'Enter degree -> '
      read (5,*,end=190,err=20) ndeg
      ndeg1=ndeg+1
      if (ndeg1.lt.1 .or. ndeg1.gt.MPOL) then
         write (0,550) 'ingrid: illegal degree'
         goto 20
      endif
      ifunpa(indexi)=ndeg
      write (0,550) 'Enter coefficients, one per prompt'
      do 675 k=0,ndeg
         write (0,1675) k
  675    read (5,*,end=190,err=671) funpar(indexr+k)
 1675 format ('  for term of power ',i1,' -> ',$)
      goto 20
*
* Execute grid generation
*
  680 continue
      write (0,551) 'Enter number of cells in x- and y-direction -> '
      read (5,*,end=190,err=20) nxc,nyc
      nxp=nxc+1
      nyp=nyc+1
      frx=1./real(nxc)
      fry=1./real(nyc)
      if (nxp.lt.2 .or. nyp.lt.2) then
         write (0,550) 'ingrid: illegal grid size'
         goto 20
      endif
      call newbuf('Enter target buffer for function',kb,
     .nxp,nyp,error,quit)
      if (quit) goto 190
      if (error) goto 20
* Prepare parameters
      b0=funpar(1)
      bx=(funpar(2)-b0)*frx
      by=(funpar(3)-b0)*fry
      bxy=(funpar(4)-b0)*frx*fry - bx*fry - by*frx
      fx=(funpar(6)-funpar(5))*frx
      fy=(funpar(8)-funpar(7))*fry
* Execute computation
      i=id(kb)
      do 682 ky=1,nyp
         y=real(ky-1)
         bili=b0+by*y
         f2=funpar(7)+fy*y
         if (ifunpa(2).ne.0) f2=sin(pi*f2)
         call poly(f2,ff2,ifunpa(4),funpar(16))
         do 682 kx=1,nxp
            x=real(kx-1)
            bili=bili+x*(bx+bxy*y)
            f1=funpar(5)+fx*x
            if (ifunpa(1).ne.0) f1=sin(pi*f1)
            call poly(f1,ff1,ifunpa(3),funpar(9))
            a(i)=bili*ff1*ff2
  682       i=i+1
      fname(kb)=' '
      proc(kb)=.false.
      used(kb)=.true.
      goto 20
*
* List function parameters
*
  690 continue
 1690 format('Corner values of bilinear: ',4f9.3)
 1691 format('Low and high values of f1: ',2f9.3)
 1692 format('Low and high values of f2: ',2f9.3)
 1693 format ('Degree of polomial for F',i1,': ',i2,'. Coefficients:')
 1694 format (2x,7f11.5)
*
      write (0,1690) (funpar(k),k=1,4)
      write (0,1691) (funpar(k),k=5,6)
      if (ifunpa(1).ne.0) write (0,550)
     .'  sin(pi * f1) will be used instead of f1'
      write (0,1692) (funpar(k),k=7,8)
      if (ifunpa(2).ne.0) write (0,550)
     .'  sin(pi * f2) will be used instead of f2'
      ndeg1=ifunpa(3)
      ndeg2=ifunpa(4)
      write (0,1693) 1,ndeg1
      write (0,1694) (funpar(k),k=9,9+ndeg1)
      write (0,1693) 2,ndeg2
      write (0,1694) (funpar(k),k=16,16+ndeg2)
      goto 600
*
* Option 8: Various utilities
*
  700 continue
      call menuh(itest,optio(1,8),quit)
      if (quit) goto 190
      goto (20,710,370,440,730,740,750),itest+1
*
* Print statistics
*
  710 continue
      call oldbuf('Enter buffer nr.',kb,.false.,error,quit)
      if (quit) goto 190
      if (error) goto 20
      zmean=0d0
      zrms=0d0
      zmin=1e15
      zmax=-1e15
      nspace=nx(kb)*ny(kb)
      id1=id(kb)-1
      l=0
      do 715 i=1,nspace
         z=a(id1+i)
         if (abs(z).lt.1e20) then
            zmin=amin1(zmin,a(id1+i))
            zmax=amax1(zmax,a(id1+i))
            zmean=zmean+a(id1+i)
            zrms=zrms+a(id1+i)**2
            l=l+1
         endif
  715 continue
      zmean=zmean/l
      zrms=dsqrt(zrms/l)
      write (0,1710) kb,zmin,zmean,zmax,zrms,l,
     .dsqrt(zrms**2-zmean**2),xmin(kb),xmax(kb),ymin(kb),ymax(kb)
 1710 format (/'Statistics for buffer nr. ',i3/6x,'Minimum',f10.3,
     .t40,'Mean',f10.3/6x,'Maximum',f10.3,t41,'RMS',f10.3/
     .'Proper points',i10,t30,'RMS about mean',f10.3/
     .'X-range ',2f8.2,t30,'Y-range ',2f8.2)
      goto 20
*
* Alternate file formats
*
  730 continue
      stdfmt=.not.stdfmt
      if (stdfmt) then
         write (0,550) 'File format fixed to default (GRID)'
      else
         write (0,550) 'Alternative file formats enabled'
      endif
      goto 20
*
* Disable or enable buffer listing
*
  740 continue
      listen=.not.listen
      if (listen) then
         write (0,550) 'Buffer listing enabled'
      else
         write (0,550) 'Buffer listing disabled'
      endif
      goto 20
*
* System commands
*
  750 continue
      write (0,550) 'Enter any number of system commands as soon as',
     .'the $ prompt appears. Exit by once hitting Ctrl-D'
      call system('sh')
      goto 20
*
* File open error
*
* 1350 call perror('ingrid')
 1350 write(0,*) "ingrid: error opening file, " // fnmenu(:l)
      goto 20
*
* Program aborts always by 'goto 999'
*
  999 end
