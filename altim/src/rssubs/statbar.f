**STATBAR -- Print status bar
*+
      SUBROUTINE STATBAR (FNC, N, TEXT)
      INTEGER*4     FNC, N
      CHARACTER*(*) TEXT
*
* This routine prints a status bar during some processing loop, keeping
* track of the progress of the processing. After each 2% of the processing
* a # is drawn, until 100% of the processing is performed.
* Upon initialisation, the maximum count to be reached should be specified.
* After that, just the continuing count should be given and the routine
* will update the status bar.
* You may specify a header in the form
*
*        2%  .    .    .    .   50%   .    .    .    .  100%
* to be printed, plus a short text (max 25 characters).
*
* Arguments:
* FNC  (input): Function:
*                0 = Initialise and print header
*               -1 = Initialise but do not print header
*               >0 = Update status bar
* N    (input): Upon initialisation: maximum count
*               During processing  : count
* TEXT (input): Text to be printed in front of status bar
*-
*  1-Dec-1993 - Created
* 21-Mar-1994 - Maximum count added in output
* 17-Oct-1994 - Proper initialisation
*-----------------------------------------------------------------------
      integer*4 n_tot,k_old,k_new,k
      save n_tot,k_old
      if (fnc.le.0) then
         if (fnc.eq.0) write (*,1000) n
	 write (*,1010) text
	 n_tot=n
	 k_old=0
	 return
      endif
      k_new=int(n*50d0/n_tot+1d-6)
      do k=k_old+1,k_new
	 write (*,1020)
      enddo
      if (n.ge.n_tot) write (*,550)
      k_old=k_new
      return

550   format(a)
1000  format(/i23,3x,
     |'2%  .    .    .    .   50%   .    .    .    .  100%')
1010  format(1a25,' ',$)
1020  format('#',$)

      end
