lgrid is a program to list the DEOS binary grid files

syntax:  lgrid [options] gridname

with optional arguments:
  fmt=fmt  : specify format for values, default=10F8.3 
  inv=inv  : specify invalid value, default= 999.999  
required arguments:
  gridname : filename of the grid

Output format (to standard output):
Header: (6F10.4,2I6) XMIN,XMAX,YMIN,YMAX,DX,DY,NX,NY
Body  : ((10F8.3 )) (GRID(K),K=1,NX*NY)
