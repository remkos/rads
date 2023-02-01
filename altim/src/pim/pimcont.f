**PIMCONT -- Doing contours in PIM
*
      subroutine pimcont(c,contours,labels)
      include "pim.inc"
      integer ncont,maxcont,i,nc
      logical contours,labels
      parameter (maxcont=1000)
      real c(*),cont(maxcont)
      external plotcont
      external lablcont

      include "pimcont.inc"
      real xtrc,ytrc

      if (.not.(contours.or.labels)) return

      xblc=xc0
      xtrc=xc1
      yblc=yc0
      ytrc=yc1
      dxp=(xtrc-xblc)/(ncx-1)
      dyp=(ytrc-yblc)/(ncy-1)

      ncont=nint(abs((rmaxc-rminc)/rintc)+1)
      if (ncont.gt.maxcont) stop 'pim: too many contours (> 1000)'
      do i=1,ncont
	 cont(i)=rminc+(i-1)*rintc
      enddo

      if (contours) then
	 write (0,550) 'Plotting contours ...'
         call pgconx(c,ncx,ncy,1,ncx,1,ncy,cont,-ncont,plotcont)
      endif

      if (labels) then
	 write (0,550) 'Plotting contour labels ...'
	 pgcmin=pgcint/2
         do i=1,ncont
	    call pgnumb(nint(cont(i)*1e3),-3,0,pgclab,nc)
	    call pgconx(c,ncx,ncy,1,ncx,1,ncy,cont(i),1,lablcont)
         enddo
      endif
550   format(a)
      end

**PLOTCONT -- External function for plotting contours
*
      subroutine plotcont(visible,xgrid,ygrid,z)
      real*4   xx,yy,z,xgrid,ygrid
      integer  visible
      include "pimcont.inc"
*
      xx=xblc+(xgrid-1)*dxp
      yy=yblc+(ygrid-1)*dyp
      call pmconv(1,xx,yy)
      if (z.ge.1e10) visible=0

      if (visible.eq.1) then
	 call pgdraw(xx,yy)
      else
         call pgmove(xx,yy)
      endif
      end

**LABLCONT -- External function for labeling contours
*
      subroutine lablcont(visible,xgrid,ygrid,z)
      real*4   z,xgrid,ygrid
      integer  visible
      include "pimcont.inc"
*
      REAL     XX, YY, XC, YC, XV1, XV2, YV1, YV2, XL, XR, YB, YT
      REAL     XN, YN
      REAL     ANGLE, XO, YO, XP, YP, DINDX, DINDY, XBOX(4), YBOX(4)
      INTEGER  I, TB
      SAVE     I
      DATA     I /0/

      xx=xblc+(xgrid-1)*dxp
      yy=yblc+(ygrid-1)*dyp
      call pmconv(1,xx,yy)
      if (z.ge.1e10) visible=0

      IF (VISIBLE.EQ.0) THEN
C        -- start of contour: reset segment counter
         I = 0
      ELSE
C        -- increment segment counter and check whether this
C           segment should be labelled
         I = MOD(I+1,PGCINT)
         IF (I.EQ.PGCMIN) THEN
C           -- find center of this segment (XC, YC)
            CALL PGQPOS(XP, YP)
            XC = (XX+XP)*0.5
            YC = (YY+YP)*0.5
C            -- find slope of this segment (ANGLE)
            CALL PGQVP(1, XV1, XV2, YV1, YV2)
            CALL PGQWIN(XL, XR, YB, YT)
            ANGLE = 0.0
            IF (XR.NE.XL .AND. YT.NE.YB) THEN
               DINDX = (XV2 - XV1) / (XR - XL)
               DINDY = (YV2 - YV1) / (YT - YB)
               IF (YY-YP.NE.0.0 .OR. XX-XP.NE.0.0)
     :           ANGLE = 57.3*ATAN2((YY-YP)*DINDY, (XX-XP)*DINDX)
            END IF
C           -- check whether point is in window
            XN = (XC-XL)/(XR-XL)
            YN = (YC-YB)/(YT-YB)
            IF (XN.GE.0.0 .AND. XN.LE.1.0 .AND.
     :          YN.GE.0.0 .AND. YN.LE.1.0) THEN
C              -- save current text background and set to erase
               CALL PGQTBG(TB)
               CALL PGSTBG(0)
C              -- find bounding box of label
               CALL PGQTXT(XC, YC, ANGLE, 0.5, PGCLAB, XBOX, YBOX)
               XO = 0.5*(XBOX(1)+XBOX(3))
               YO = 0.5*(YBOX(1)+YBOX(3))
C              -- plot label with bounding box centered at (XC, YC)
               CALL PGPTXT(2.0*XC-XO, 2.0*YC-YO, ANGLE, 0.5, PGCLAB)
C              -- restore text background
               CALL PGSTBG(TB)
            END IF
         END IF
      END IF
      CALL PGMOVE(XX,YY)
      END
