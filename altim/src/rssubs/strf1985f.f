**STRF1985F -- Construct date string from seconds since 1985
*+
      SUBROUTINE STRF1985F (STRING, SEC)
      CHARACTER(len=26) STRING
      REAL*8 SEC

* This routine formats relative time in seconds since 1985 into a
* character string of 26 bytes in the form YYYY-MM-DD HH:MM:SS.SSSSSS.
*-
* Copyright (c) Remko Scharroo, Altimetrics LLC
*-----------------------------------------------------------------------
      integer isec
      isec = int(sec)
      call strf1985(string,"%Y-%m-%d %H:%M:%S",isec)
      write (string(20:26),'(f7.6)') sec-isec
      end
