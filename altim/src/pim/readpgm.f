**READPGM -- Read PGM file, save to pixel map
*+
      subroutine readpgm(filenm,buf,pixmap,nx,ny,ci)
      character*(*) filenm
      integer*4 pixmap(*),nx,ny,ci(0:255)
      integer*1 buf(*)
      character*256 dir
      integer l
*
* Read a PGM file from disk. Save the contents in a pixel map buffer,
* while re-assigning the color indices
*
* Arguments:
*  FILENM  (input): PGM filename
*  BUF    (output): Temporary working space
*  PIXMAP (output): Array of NX by NY containing assigned color indices
*  NX, NY (output): Size of the pixel map
*  CI(I)   (input): Colour index to be assigned to index I in the PGM file
*-
      integer GRIFIL,i
      character*80 line

      i=GRIFIL(filenm)
      if (i.le.0) then
         dir='/user/altim'
	 call checkenv('ALTIM',dir,l)
	 dir(l+1:)='/pim/'//filenm
         i=GRIFIL(dir)
         if (i.le.0) stop 'PGM file not found'
      endif
      call readpgm1(line,i)
      if (line(1:2).ne.'P5') stop 'Illegal PGM file'
1     call readpgm1(line,i)
      if (line(1:1).eq.'#') goto 1
      read (line,*) nx,ny
      call readpgm1(line,i)
      call GRRFIL(i,nx*ny,buf)
      call GRCFIL(i)

      call readpgm2(buf,pixmap,nx,ny,ci)
      end
