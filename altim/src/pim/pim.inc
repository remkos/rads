**file: pim.inc
*+
* Standard definitions
*
      integer*2 minint2,maxint2
      parameter (maxint2=32767,minint2=-maxint2-1)
      integer*4 minint4,maxint4
      parameter (maxint4=2147483647,minint4=-maxint4-1)

* MAXCOL = maximum number of colour indices
* MAXGRD = maximum number of grid points
*
      INTEGER MAXCOL
      PARAMETER (MAXCOL=255)
      integer maxgrd
      parameter (maxgrd=1280*1024)
      real*4 work1(maxgrd),work2(maxgrd),work3(maxgrd),work4(maxgrd),
     .	work5(maxgrd)
      COMMON /WORK/ WORK1, WORK2, WORK3, WORK4, WORK5
*
* XW0, XW1, YW0, YW1 = plot area boundaries
* XV0, XV1, YV0, YV1 = viewport settings
* PX, PY = Pixels
* PROJECT = projection type
* INTER = interpolation degree; 0=none, 1=some, 2=highest resolution
*
      REAL*4  XW0,XW1,YW0,YW1,XV0,XV1,YV0,YV1
      INTEGER PX,PY
      COMMON /SIZES/ XW0,XW1,YW0,YW1,XV0,XV1,YV0,YV1,PX,PY
      INTEGER PROJECT, INTER
      COMMON /MAP/ PROJECT, INTER
*
* NC1 = Number of colours - 1
* NC2 = Number of shades per colour - 1
* C_1 = Highest colour index
* SH_MAT = Conversion matrix INDEX,LIGHTNESS => NEW INDEX
*
      INTEGER NC1,NC2,C_1,SH_MAT(0:MAXCOL,0:MAXCOL)
      REAL RGB(3,0:MAXCOL),SHADE0,SHADE1
      COMMON /COLORS/ NC1,NC2,C_1,SHADE0,SHADE1,SH_MAT,RGB
*
* XS0, XS1, YS0, YS1 = Surface grid boundaries
* RMINS, RMAXS = Range for surface grid
* NSX, NSY = Surface grid size
* SNAME = Grid name
* SARGUM = Command line arguments
* etc. (C=Contour grid, L=Illumination grid)
*
      REAL*4  XS0,XS1,YS0,YS1,RMINS,RMAXS,RMEANS,RSIGS,
     |        XSCOG(4000),YSCOG(4000)
      INTEGER NSX,NSY
      LOGICAL CELLS,SLOG
      CHARACTER*1024 SNAME,SARGUM
      COMMON  /CGRIDS/ XS0,XS1,YS0,YS1,RMINS,RMAXS,
     |        RMEANS,RSIGS,XSCOG,YSCOG,NSX,NSY,SNAME,SARGUM,CELLS,SLOG

      REAL*4  XC0,XC1,YC0,YC1,RMINC,RMAXC,RINTC,RMEANC,RSIGC
      INTEGER NCX,NCY
      LOGICAL CELLC
      CHARACTER*1024 CNAME,CARGUM
      COMMON  /CGRIDC/ XC0,XC1,YC0,YC1,RMINC,RMAXC,
     |        RMEANC,RSIGC,NCX,NCY,CNAME,CARGUM,CELLC,RINTC

      REAL*4  XL0,XL1,YL0,YL1,RMINL,RMAXL,RMEANL,RSIGL
      INTEGER NLX,NLY
      LOGICAL CELLL
      CHARACTER*1024 LNAME,LARGUM
      COMMON  /CGRIDL/ XL0,XL1,YL0,YL1,RMINL,RMAXL,
     |        RMEANL,RSIGL,NLX,NLY,LNAME,LARGUM,CELLL

      REAL*4  RMINB,RMAXB
      COMMON  /CGRIDB/ RMINB,RMAXB
*
* POSTSCRIPT = .true. if PostScript device
* DITHER1, DITHER2 = .true. if dithering is applied in colour or shading
*   direction
* DUMPED = .true. if pixel map is dumped or not used.
* DEBUG gives extra output if .true.
*
      LOGICAL POSTSCRIPT,DITHER1,DITHER2,DUMPED,DEBUG
      COMMON /MODE/ POSTSCRIPT,DITHER1,DITHER2,DUMPED,DEBUG
*
* Set special colour indices. Colour range starts at c_0
*
      INTEGER C_BAD, C_COAST, C_LAND, C_FG, C_BG, C_CONT,
     .  C_BLUE, C_RED, C_GREY, C_WHITE, C_BLACK,
     .  C_HALF, C_LGREY, C_DGREY, C_0
      PARAMETER (
     .  C_BG=0, C_FG=1, C_LAND=2, C_BAD=3, C_COAST=4, C_CONT=5,
     .  C_WHITE=12, C_LGREY=9, C_GREY=8, C_DGREY=7,
     .  C_BLUE=10, C_RED=11, C_BLACK=6, C_HALF=15, C_0=16)
*
* Set default character height
*
      REAL DEF_CH
      PARAMETER (DEF_CH=0.8)
*
* Command line options
*
      CHARACTER*1024 NAME,ARGUM
      REAL FS,CI,LW,LS,CH,DUM,X,Y,FACT,RMIN,RMAX
      INTEGER IOS,IX,IY
      LOGICAL POPCMD, POP1CMD,POPNEXTCMD
      COMMON /CMD/ NAME,ARGUM,FS,CI,LW,LS,CH,DUM,X,Y,FACT,
     .  IOS,IX,IY,RMIN,RMAX
*
* To show trends in stacks of grids or XGF
*
      LOGICAL GRIDTRND, XGFTRND
      INTEGER NSTACK, MSTACK
      PARAMETER (MSTACK=50)
      CHARACTER*80 XNAME
      REAL    TX(MSTACK),TY(MSTACK)
      COMMON /TRENDS/ XNAME,GRIDTRND,XGFTRND,NSTACK,TX,TY
*
* Other settings
*
      CHARACTER*80 DEVICE
      INTEGER DEVICE_ID
      COMMON /SETTINGS/ DEVICE, DEVICE_ID
*-
