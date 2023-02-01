**ALTBIAS -- Altimeter range biases and drifts
*+
      SUBROUTINE ALTBIAS (PROD, UTC, COR0, COR1, COR2)
      implicit none
      INTEGER  PROD, UTC
      REAL*8   COR0, COR1, COR2
*
* This function returns altimeter range bias and drifts for various altimeter
* products. These range deltas can be either constant (range bias)
* or time variant (drift). Because of the different problems associated
* to each satellite product, the actual source of the biases and drifts
* are different.
*
* Input to the function are the satellite product identifier (PROD),
* and the epoch in UTC seconds of 1985 (UTC). The function returns
* three deltas: a constant bias (COR0) and two time-variant terms
* (COR1 and COR2), all in metres. To correct the range, COR0 should
* be SUBTRACTED FROM the ranges provided on the data and COR1 and COR2
* should be ADDED (due to conventions).
*
* Thus:    RANGE (corrected) = RANGE - COR0 + COR1 + COR2
* Or:      SSH   (corrected) = SSH   + COR0 - COR1 - COR2
*
* Arguments (see below for further description):
*  PROD      (input) : Satellite ID or Product identifier
*  UTC       (input) : Time in UTC seconds since 1.0 Jan 1985
*                      or TOPEX/POSEIDON Cycle number
*  COR0     (output) : Range bias, constant part
*  COR1     (output) : Range bias, variable part
*  COR2     (output) : USO correction
*
* Satellite ID or Product identifier (PROD)
* -----------------------------------------
* 3/300     : Geosat GDR
* 4/400     : ERS-1 OPR Version 5 and higher
* 5/500     : TOPEX (T/P Merged GDR, new format)
* 6/600     : POSEIDON (T/P Merged GDR, new format)
* 7/700     : ERS-2 OPR Version 5 and higher
* 8         : GFO GDR
* 9         : Jason-1 (I)GDR
* 10        : Envisat (I)GDR
*   401     : ERS-1   OPR Version 3
*   410/710 : ERS-1/2 QLOPR
*   411/711 : ERS-1/2 URA
*   412/712 : ERS-1/2 IGDR/RGDR
*   420/720 : ERS-1/2 WAP
*   501/601 : TOPEX/POSEIDON Merged GDR, old format
*
* Meaning of COR0, COR1, and COR2
* -------------------------------
*           COR0            COR1                   COR2
* Geosat    range bias    Internal calibration     USO drift
* ERS-1     range bias    SPTR bias                USO drift
* TOPEX     range bias    Internal calibration     USO drift correction
* POSEIDON  range bias         0                       0
* ERS-2     range bias    SPTR bias                USO drift
* GFO            0             0                       0
* Jason-1        0             0                       0
* Envisat   range bias         0                   USO drift
*
* References
* - The GEOSAT range bias is taken to be 12.4 cm (Brian Beckley, priv.
*   comm., 2002)
* - The ERS-1 range bias was determined during the Venice Arc Calibration
*   to be -41.6 cm (measuring short).
*   C. R. Francis et al., "The Calibration of the ERS-1 Radar Altimeter --
*   The Venice Calibration Campaign", ESA Report ER-RP-ESA-RA-0257, issue
*   2.0, 1 March 1993.
*   Nevertheless, we adopt the value of -40.92 cm as measured on ERS-2.
*   (see also Stum et al. [1998], below)
* - The ERS-2 range bias is determined prelaunch at Dornier to be -40.92 cm
*   (E. Schied). This value is included in the OPR product, so should lead
*   to a zero range bias. This correction is not applied in the FD data
*   (URA/IGDR/RGDR/URA) and is thus provided off-line.
* - The impact of the change from OPR version 3 to OPR version 6 results
*   in a change of range of 24 mm (to be subtracted from the range)
*   J. Stum, F. Ogor, P-Y. Le Traon, J. Dorandeu, P. Gaspar, J-P. Dumont,
*   "An inter-calibration study of TOPEX/POSEIDON, ERS-1 and ERS-2
*   altimetric missions", Final Report of IFREMER contract No. 97/2 426 086/C,
*   Ref. CLS/DOS/NT/98.070, CLS, Ramonville, France, April 1998
* - The ERS-1/2 SPTR bias jumps were recovered by ESTEC (M. Roca) and
*   tabled for each period inbetween altimeter anomalies.
*   C. R. Francis and M. Roca, "Spontaneous Jumps in the Relative Bias of
*   ERS-1 and ERS-2 Altimetry: Their Origin and Cure", AGU 96 Spring
*   Meeting, Baltimore, 20-24 May 1996.
* - New SPTR and USO tables were introduced late 2000. They are the
*   default since 14 Feb 2001. USO values remained the same, but the
*   SPTR values are based on a new retrieval algorithm:
*   A. Martini and P. Femenias, "The ERS SPTR2000 Altimeric Range
*   Correction: Results and Validation", ESA ERE-TN-ADQ-GSO-6001, Nov 2000.
* - The ERS-1/2 USO range correction is due to changes in the Oscilator
*   frequency which are not accounted for in the ERS-1 altimeter products.
*   Note that for OPR and QLOPR/URA/IGDR different reference frequencies
*   have been used, so the USO range drift correction is different for
*   these products.
*   C. Loial and J. Benveniste, ESRIN, priv. comm.
* - The ERS-1/2 Internal calibration correction is always applied to the
*   data.
* - The TOPEX Internal calibration is registered at Wallops, and consists
*   of actual changes in instrument and temperature effects. Here,
*   we give the combined effect.
*   G. S. Hayne, D. W. Hancock III, C. L. Purdy, "TOPEX Altimeter Range
*   Stability Estimates from Calibration Mode Data", TOPEX/POSEIDON
*   Research News, Issue 3, pp. 18-20, October 1994.
*   G. S Hayne, "TOPEX Altimeter Range Stability Estimate Update",
*   http://topex.wff.nasa.gov/
* - The TOPEX USO drift was applied to the old GDR data with the wrong sign.
*   Therefore, an external table is available to correct for this error.
*   D. W. Hancock III and G. S Hayne, "Error in TOPEX Oscillator Drift
* - After correction the TOPEX range bias should effectively be zero.
*   The same is assumed for POSEIDON. This is done so for the new format.
* - Correct Envisat RA-2 for USO drift according to a fit to the ESA USO
*   drift tables. Ref: "USO Clock Period Corrections",
*   http://earth.esa.int/pcs/envisat/ra2/auxdata/
*   Also add the same 409.2 mm bias as ERS-1/2.
*-
* $Log: altbias.f,v $
* Revision 1.9  2011/05/12 21:25:18  rads
* - Skip commented out lines in ERS-1/2 SPTR files
*
* Revision 1.8  2007/02/21 18:26:30  rads
* - Small changes to work better with gfortran (gcc 4.2)
*
* Revision 1.7  2006/01/20 02:54:15  rads
* - Final version of Envisat USO drift function. Use only for CMA < 7.1
*
* Revision 1.5  2005/08/01 14:34:50  rads
* - New function for Envisat USO correction
* - ERS-2 WAP biases defined
*
* Revision 1.4  2005/01/19 03:13:49  remko
* - Renewed Envisat USO drift function based on actual data
*
* Revision 1.3  2004/09/22 18:55:17  remko
* - Limit N1 RA-2 USO drift to dates AFTER 13 Jun 2003
*
* Revision 1.2  2004/09/21 16:15:26  remko
* - Removed ERSDRIFT. Only ERSDRIFT2K currently used.
* - Added Envisat USO drift and bias
* - Included dummies (0 bias and drift) for other altimeters
*
* 19-Jun-2002 - Patched error in ERS-1 OPR OSU corrections. Affects
*               period 13-Oct-1991 to 21-Mar-1995.
* 28-Feb-2002 - Allow linear extrapolation of GEOSAT at end of table
* 25-Feb-2002 - Implemented linear interpolation into GEOSAT internal
*               cal and USO tables; 12.4 cm (i.s.o. 17.5 cm) range bias
* 21-Dec-2000 - New SPTR and USO tables for ERS-1 and -2
* 10-Dec-1998 - Implemented extra 24 mm bias for ERS-1 OPR v3
*  9-Dec-1998 - Check on premature end (i.e. no blank line) when
*               reading USO or SPTR table
*  4-Dec-1998 - Ignore lines with "****" in USO table
* 11-Nov-1998 - Increased table end by 7 days
* 23-Oct-1997 - Removed check on SPTR table end
* 11-Jul-1997 - Added TOPEX and POSEIDON new format (bias=0)
* 11-Jun-1997 - Added common block in ERSDRIFT for exporting
*	        One-too-high index in ERSDRIFT removed
* 28-Feb-1997 - Bug caused USO correction to be zero for Cycle 132
* 30-Oct-1996 - Sign conventions changed. ERS-1 bias changed.
* 18-Oct-1996 - New SPTR correction format
* 24-Jul-1996 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4 satel,iprod
      real*8	t,t0,t1,t2,t3,a0,a1,b0,b1

      iprod=abs(prod)
      if (iprod.lt.100) then
	 satel=iprod
	 iprod=satel*100
      else
         satel=iprod/100
      endif

      if (satel.eq.3) then

* GEOSAT: constant bias of 12.4 cm; HCAL and USO from table

	 cor0=12.4d-2
	 call geodrift(utc,cor1,cor2)

      else if (satel.eq.4) then

* ERS-1: range bias is -40.92 cm - SPTR - USO drift
* For URA data correct USO drift by +0.200 Hz => -10.60 mm
* For WAP data correct USO drift by +0.150 Hz =>  -7.95 mm
* For OPR version 3 data add 24 mm to the bias [Stum et al., 1998]

	 cor0=-40.92d-2
	 call ersdrift2k(1,utc,cor1,cor2)
	 if (iprod  .eq.401) cor0=cor0+24d-3
	 if (iprod/10.eq.41) cor2=cor2-10.60d-3
	 if (iprod/10.eq.42) cor2=cor2-7.95d-3

      else if (iprod.eq.501) then

* TOPEX: range bias is 0.0 cm - Internal Cal - USO drift error

	 cor0=0d0
	 call tpxdrift(utc,cor1,cor2)

      else if (satel.eq.7) then

* ERS-2 OPR: range bias is   0.0  cm - SPTR - USO drift
*      RGDR: range bias is -40.92 cm - SPTR - USO drift
* For URA data correct USO drift by -0.040 Hz => +2.12 mm
* For WAP data correct USO drift by -0.090 Hz => +4.77 mm

	 call ersdrift2k(2,utc,cor1,cor2)
	 cor0=0d0
	 if (iprod  .eq.712)  cor0=-40.92d-2
	 if (iprod/10.eq.71)  cor2=cor2+2.12d-3
	 if (iprod/10.eq.72) then
	    cor0=234d-3
	    cor2=cor2+4.77d-3
	 endif

      else if (satel.eq.10) then

* Envisat: apply same range bias as ERS-1 plus USO drift

	 t=(utc/86400-5478)/365.25d0	! years since 2000.0
	 a0 =32.7587d0		! Values determined up to Cycle 41
	 a1 =-2.1288d0		! from auxiliary USO files
	 b0 =-3.45023d0
	 b1 =-3.73254d0
	 t0 = 3.48070d0
	 t1 = 4.31319d0
	 t2 = 4.7113d0
	 t3 = 5.5051d0
	 cor0=-40.92d-2
	 cor1=0d0
	 if (t.lt.t0) then
	    cor2=a0
	 else if (t.lt.t1) then
	    cor2=a0+b0*(t-t0)
	 else if (t.lt.t2) then
	    cor2=a0+b0*(t-t0)+b1*(t-t1)
	 else if (t.lt.t3) then
	    cor2=a0+b0*(t-t0)+b1*(t2-t1)
	 else
	    cor2=a0+a1+b0*(t-t0)
	 endif
	 cor2=-1d-3*cor2	! Change sign and go to meters

      else

* All others: no biases or drifts. This includes:
* - POSEIDON and TOPEX new format
* - GFO
* - Jason-1

	 cor0=0d0
	 cor1=0d0
	 cor2=0d0
      endif
      return

      end

**ERSDRIFT2K -- Determine SPTR and USO drift for ERS-1 and -2 (v.2000)
*+
      SUBROUTINE ERSDRIFT2K (ERS, UTC, SPTR, USO)
      INTEGER*4 ERS, UTC
      REAL*8    SPTR, USO
*
* This routine extracts the SPTR correction and USO drift correction
* to the ERS-1 and ERS-2 ranges from off-line tables (version 2000).
* The SPTR correction and the USO drift correction
* both have to be ADDED to the altimeter range.
* Both SPTR and USO corrections are given in metres and apply
* to URA data. Additional corrections are to be applied for OPR data
* (ERS-1: add 10.60 mm; ERS-2: subtract 2.12 mm)
*
* SPTR:
* Values change stepwise after each time the altimeter was switched off.
* When an epoch after the last table entry is requested, the last value is
* returned. When no number is available, i.e. when no SPTR
* measurement was made between two anomalies, the average SPTR value is
* returned. This is +20 mm for ERS-1 and -20 mm for ERS-2.
* The SPTR2K values are more precise than the older values
* obtained by ERSDRIFT.
*
* OSU:
* The new USO tables (version 2000) are the same as the older ones.
* However, in ERSDRIFT2K linear inter/extrapolation in stead of stepwise
* change is implemented. The differences are marginal.
*
* On the first call the tables are loaded. On subsequent calls, the
* values are read from memory.
*
* Arguments:
* ERS   (input) : 1=ERS-1, 2=ERS-2
* UTC   (input) : Time in UTC seconds of 1985
* SPTR (output) : SPTR correction (metres)
* USO  (output) : USO drift correction (metres)
*-
* 21-Dec-2000 - Created by Remko Scharroo
*-----------------------------------------------------------------------
      integer*4 max_index,t_max_index(4)
      parameter (max_index=3000)
      integer*4 t_index(4),t_utc(max_index,4),t0(4)
      real*8    t_rng(max_index,4)
      logical   first/.true./
      integer   j,k,l,unit,freeunit,yymmdd,hh,mm,ss,day,nr,flag
      real*8    delta,sec85,sig
      character*120 dirname,filename
      character*1 dum

      common /cersdrift2k/ t_utc,t_rng,t_max_index

      save

600   format (a,'/data/bias/ers',i1,'_ra_',a,'.txt')
610   format (a1,i6,3(1x,i2),2(f10.1),i11,i6)
620   format (11x,3(1x,i2),4x,i5,35x,f8.3)
      if (first) then
	 first=.false.

* Load set date offsets

	 t0(1)=157766400	! 1.0 Jan 1990 - 1.0 Jan 1985
	 t0(2)=t0(1)
	 t0(3)=(48454-46066)*86400	! Launch of ERS-1
	 t0(4)=(49828-46066)*86400	! and ERS-2

* Read ERS SPTR correction files

	 unit=freeunit()
	 dirname='/user/altim'
	 call checkenv('ALTIM',dirname,l)
	 do j=1,2
	    write (filename,600) dirname(:l),j,'sptr'
            open (unit,file=filename,status='old')
	    k=0
30	    read (unit,610,end=39) dum,yymmdd,hh,mm,ss,delta,sig,nr,flag
	    if (dum.eq.'#') goto 30 ! Skip commented out lines
	    k=k+1
	    if (k.gt.max_index)
     |		call fin("ERSDRIFT2K: too many SPTR table entries")
	    if (flag.eq.0) then
	    ! Use delta value if flag=0, otherwise use:
	    else if (j.eq.1) then
	       delta=+20	! Average value for ERS-1
	    else
	       delta=-20	! Average value for ERS-2
	    endif
	    t_utc(k,j)=sec85(4,yymmdd*1d6+hh*1d4+mm*1d2+ss)
	    t_rng(k,j)=delta/1d3
*	    write (0,*) 'sptr:',k,t_utc(k,j),t_rng(k,j)
	    goto 30
39	    continue
	    t_utc(1,j)=-2147483647-1
	    t_utc(k+1,j)=2147483647
	    t_max_index(j)=k
	    t_index(j)=(k+1)/2
	 enddo
	 close (unit)

* Read ERS USO correction files

	 do j=3,4
	    write (filename,600) dirname(:l),j-2,'uso'
            open (unit,file=filename,status='old')
	    k=1
20	    read (unit,620,end=29) hh,mm,ss,day,delta
	    k=k+1
	    if (k.gt.max_index)
     |		call fin("ERSDRIFT2K: too many USO table entries")
	    if (hh.eq.99) then
* When hh:mm:ss = 99:99:99, use the average epoch of 11:00:00
	       hh=11; mm=00; ss=00
	    endif
	    t_utc(k,j)=t0(j)+day*86400+hh*3600+mm*60+ss
	    t_rng(k,j)=delta/1d3
*	    write (0,*)	'uso:',k,t_utc(k,j),t_rng(k,j)
	    goto 20
29	    continue
	    t_utc(1,j)=-2147483647-1
	    t_utc(k+1,j)=2147483647
	    t_max_index(j)=k
	    t_index(j)=(k+1)/2
	 enddo
	 close (unit)

      endif

* Scan for correct period in SPTR and USO table

      do j=ers,4,2
110	 if (utc.ge.t_utc(t_index(j)+1,j)) then
	    t_index(j)=t_index(j)+1
	    goto 110
	 else if (utc.lt.t_utc(t_index(j),j)) then
	    t_index(j)=t_index(j)-1
	    goto 110
	 endif
      enddo

* SPTR: take previous value

      SPTR=t_rng(t_index(ers),ers)

* USO: linear inter/extrapolation between previous and next value

      j=ers+2
      k=min(t_index(j),t_max_index(j)-1)
      USO=t_rng(k,j)+(t_rng(k+1,j)-t_rng(k,j))*(utc-t_utc(k,j))
     |		/(t_utc(k+1,j)-t_utc(k,j))

      end

**TPXDRIFT -- Determine Internal Cal and USO drift for TOPEX
*+
      SUBROUTINE TPXDRIFT (UTC, CAL, USO)
      INTEGER*4 UTC
      REAL*8    CAL, USO
*
* This routine extracts the internal calibration correction (Wallops
* correction) and USO drift error correction to the TOPEX ranges from
* off-line tables. The internal calibration correction and the USO
* drift correction are to be ADDED to the altimeter range.
* Both CAL and USO corrections are given in metres.
*
* On the first call the tables are loaded. On subsequent calls, the
* values are read from memory.
*
* Arguments:
* UTC  (input) : Time in UTC seconds of 1985
* CAL (output) : Internal calibration correction (metres)
* USO (output) : USO drift error correction (metres)
*-
* 24-Jul-1996 - Created by Remko Scharroo
* 25-Jul-1996 - T0 was wrong by half a day.
*-----------------------------------------------------------------------
      integer*4 max_index,u_max_index,s_max_index
      parameter (max_index=999)
      real*8    u_rng(max_index),s_rng(max_index)
      logical   first/.true./
      integer   cyc,i,l,unit,freeunit
      real*8    delta,t0,topex_length
      parameter (topex_length=856711.542d0,t0=242977173d0)
      character*80 line,dirname

      save

      if (first) then
	 first=.false.

* Read TOPEX Internal Calibration Table

	 unit=freeunit()
	 dirname='/user/altim'
	 call checkenv('ALTIM',dirname,l)
         open (unit,file=dirname(:l)//'/data/bias/RngStbUp.txt',
     |          status='old')
	 call headscan(unit,'Cyc  Count',0)
50	 read (unit,550,end=59) line
	 if (line.eq.' ' .or. line(:1).eq.'<') goto 59
	 read (line,54) cyc,delta
54	 format (i3,9x,f7.3)
	 s_rng(cyc)=-delta/1d3
*	 write (0,*) 'cal:',cyc,s_rng(cyc)
	 goto 50
59	 continue
	 s_max_index=cyc
	 close (unit)

* Read TOPEX Oscillator Drift Table

         open (unit,file=dirname(:l)//'/data/bias/OscDrift.txt',
     |          status='old')
	 call headscan(unit,'Cycle midpoint',2)
40	 read (unit,550,end=49) line
	 if (line.eq.' ') goto 49
	 if (line(64:68).eq.'SSALT') then
	    read (line,44) cyc
	    delta=1d30
	 else
	    read (line,44) cyc,delta
         endif
44	 format (i5,57x,f6.2)
	 u_rng(cyc)=delta/1d3
*	 write (0,*) 'uso:',cyc,u_rng(cyc)
	 goto 40
49	 continue
	 do i=cyc+1,max_index
	    u_rng(i)=0
	 enddo
	 u_max_index=max_index
	 close (unit)

      endif

* Convert to cycle number, then pick number from table

      if (utc.lt.1000) then
	 cyc=utc
      else
	 cyc=int((utc-t0)/topex_length)
      endif

      if (cyc.lt.1)
     |  call fin("TPXDRIFT: Requested epoch before table range")
      if (cyc.gt.u_max_index .or. cyc.gt.s_max_index)
     |	call fin("TPXDRIFT: Requested epoch beyond table range")
      if (u_rng(cyc).gt.1d20)
     |	write (0,550) "TPXDRIFT: Requested epoch is not TOPEX"
      
      CAL = s_rng(cyc)
      USO = u_rng(cyc)

550   format (a)
      end

**GEODRIFT -- Return GEOSAT calibration and USO drift values

      SUBROUTINE GEODRIFT (UTC, CAL, USO)
      INTEGER*4 UTC
      REAL*8    CAL, USO
*
* This routine extracts the internal calibration correction correction
* and USO drift error correction (Hayne and Hancock) from an off-line
* table. The internal calibration correction and the USO
* drift correction are to be ADDED to the altimeter range.
* Both CAL and USO corrections are given in metres.
*
* On the first call the table is loaded. On subsequent calls, the
* values are read from memory.
*
* Arguments:
* UTC  (input) : Time in UTC seconds of 1985
* CAL (output) : Internal calibration correction (metres)
* USO (output) : USO drift error correction (metres)
*-
* 24-Oct-1997 - Created by Remko Scharroo
* 20-Dec-2000 - Added SAVE command
* 21-Feb-2002 - Implemented linear interpolation
* 28-Feb-2002 - Allow linear extrapolation at end of table
*-----------------------------------------------------------------------
      integer*4 n,nmax,unit,freeunit,sec
      parameter (nmax=1674)
      real*8 value(2,nmax),x
      logical first /.true./
      character*80 dirname
      integer*4 l

      save

      if (first) then
         first=.false.
         unit=freeunit()
	 dirname='/user/altim'
	 call checkenv('ALTIM',dirname,l)
	 open (unit,file=dirname(:l)//'/data/bias/hcal_uso.txt',
     |		status='old')
	 do n=1,3
	    read (unit,*)
	 enddo
	 do n=1,nmax
    	    read (unit,*) sec,cal,uso
	    value(1,n)=cal/1d3
	    value(2,n)=uso/1d3
         enddo
         close (unit)
      endif

      x=utc/86400d0-76d0
      n=min(int(x),nmax-1)
      x=x-n
      if (n.lt.1 .or. n.gt.nmax-1)
     |	call fin("GEODRIFT: Time outside table")
      cal=value(1,n)*(1-x)+value(1,n+1)*x
      uso=value(2,n)*(1-x)+value(2,n+1)*x
      
      return

      end

      subroutine headscan(unit,text,skip)
      integer unit,skip,i
      character text*(*),line*132

10    read (unit,550,err=1300,end=1300) line
      if (index(line,text).ne.0) then
         do i=1,skip
            read (unit,*)
         enddo
         return
      else
         goto 10
      endif

550   format (a)

1300  call fin("HEADSCAN: Error scanning for table header")
      end
