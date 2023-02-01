**PMQINF -- return general PMPLOT information
*+
      SUBROUTINE PMQINF (NAME, VALUE)
      CHARACTER*(*) NAME
      REAL VALUE(*)
*
* Query some of the PMPLOT settings such as version-number, window settings,
* etc. NAME is a string identifying which information should be returned.
* This character value can be given both in lower and upper case.
* VALUE is a buffer of returns values.
*
* Arguments:
*  NAME   (input) : Type of information to query (case is irrelevant):
*       'VERSION' : Version number will be returned.
*          'AREA' : Min/Max longitude and Min/Max latitude in plot.
*      'USERAREA' : Area set by the user in PMDEF.
*         'SCALE' : Scale (1:SCALE) of the plot.
*  VALUE (output) : Value(s) returned.
*--
*  6-Apr-1993 - Created [Remko Scharroo]
*-----------------------------------------------------------------------
      INCLUDE 'pmplot.inc'
      CHARACTER*32 CHR
      CALL GRTOUP(CHR,NAME)
      IF (NAME.EQ.'VERSION') THEN
         VALUE(1)=PVERSION
      ELSE IF (NAME.EQ.'SCALE') THEN
         VALUE(1)=PSCALE
      ELSE IF (NAME.EQ.'AREA') THEN
         VALUE(1)=LONMIN
         VALUE(2)=LONMAX
         VALUE(3)=LATMIN
         VALUE(4)=LATMAX
      ELSE IF (NAME.EQ.'USERAREA') THEN
         VALUE(1)=LONUMIN
         VALUE(2)=LONUMAX
         VALUE(3)=LATUMIN
         VALUE(4)=LATUMAX
      ENDIF
      END
