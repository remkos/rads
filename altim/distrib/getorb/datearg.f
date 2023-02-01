**DATEARG -- Processes standard datation arguments
*+
      FUNCTION DATEARG (ARG, T0, T1, DT)
      LOGICAL*4 DATEARG
      CHARACTER*(*) ARG
      REAL*8    T0, T1, DT

* This function processes standard datation arguments of either of the
* following forms:
*
*   mjd=t0[,t1[,dt]] : Modified Julian Dates
*   d50=t0[,t1[,dt]] : Days since 1.0 Jan 1950
*   d00=t0[,t1[,dt]] : Days since 1.0 Jan 2000
*   sec=t0[,t1[,dt]] : Seconds since 1.0 Jan 1985
*   s85=t0[,t1[,dt]] : Seconds since 1.0 Jan 1985
*   s00=t0[,t1[,dt]] : Seconds since 1.0 Jan 2000
*   ymd=t0[,t1[,dt]] : [YY]YYMMDD.DDD or [YY]YYMMDDHHMMSS.SSS
*   doy=t0[,t1[,dt]] : YYDDD.DDD
*     t=t0[,t1[,dt]] : MJD.DDD or [YY]YYMMDD.DDD or [YY]YYMMDDHHMMSS.SSS
*
* Normally ARG is read from the argument line of a command using the GETARG
* routine. When ARG is parsed through to DATEARG it checks whether ARG is
* of one of the above forms. If not, DATEARG gets the value .FALSE.
* If ARG is of one of the standard datation forms, the values of T0 and
* T1 (when given) are read and interpreted as stated above and converted
* to Seconds since 1.0 Jan 1985. When provided, DT is read too, but not
* converted. DATEARG will get the value .TRUE.
*
* When T1 or DT are not provided in ARG, the arguments of the function call
* will keep their values unchanged.
*
* Arguments:
*   DATEARG (output): .TRUE. if ARG is in a standard datation format
*   ARG      (input): datation argument (see above)
*   T0, T1  (output): datation arguments converted to SEC85
*                     (Seconds since 1.0 Jan 1985)
*   DT      (output): Third datation argument (not converted)
*-
* 03-Nov-1999 : Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4 ios,l,mode
      real*8 sec85,tt0,tt1,dtt

* Initialize the values to an unlikely value

      tt0 = 1d50
      tt1 = 1d50
      dtt = 0d0

* Scan the argument for datation format

      l=index(arg,'=')
      if (l.ne.2 .and. l.ne.4) then
         datearg = .false.
	 return
      else if (arg(:l).eq.'t=' .or. arg(:l).eq.'T=') then
	 mode = 0
      else if (arg(:l).eq.'mjd=' .or. arg(:l).eq.'MJD=') then
         mode = 1
      else if (arg(:l).eq.'doy=' .or. arg(:l).eq.'DOY=') then
         mode = 3
      else if (arg(:l).eq.'ymd=' .or. arg(:l).eq.'YMD=') then
         mode = 5
      else if (arg(:l).eq.'sec=' .or. arg(:l).eq.'SEC=' .or.
     |         arg(:l).eq.'s85=' .or. arg(:l).eq.'S85=') then
         mode = -1
      else if (arg(:l).eq.'s00=' .or. arg(:l).eq.'S00=') then
         dtt = 473299200d0
         mode = -1
      else if (arg(:l).eq.'d00=' .or. arg(:l).eq.'D00=') then
         dtt = 51544d0
	 mode = 1
      else if (arg(:l).eq.'d50=' .or. arg(:l).eq.'D50=') then
         dtt = 33282d0
	 mode = 1
      else
         datearg=.false.
	 return
      endif

* Read the values, and convert when specified

      read (arg(l+1:),*,iostat=ios) tt0,tt1,dt
      if (abs(tt0).lt.1d49) t0 = sec85(mode,tt0+dtt)
      if (abs(tt1).lt.1d49) t1 = sec85(mode,tt1+dtt)
      datearg=.true.
      end
