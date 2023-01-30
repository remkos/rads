**STRIP -- Strip the command A and argument B from string A
*+
      SUBROUTINE STRIP (A, B)
      CHARACTER*(*) A, B
*
* This routine takes string A and strips of the first word. This is
* then stored in B while the remaining string is stored back in B.
* If A is empty, A and B will be returned unaltered.
* If A contains one string, B will be A and A will be empty.
*
* Arguments:
*  A  (input): String of words
*    (output): String of words minus the first word
*  B (output): First word in A
*
* Example:
* A = 'The big brown fox' => STRIP => B = 'The', A = 'big brown fox'
*-
      integer i,l
      l=len(a)

10    i=index(a,' ')
      if (a.eq.' ') then
	 return
      else if (i.eq.0 .or. i.eq.l) then
	 b=a
	 a=' '
      else if (i.eq.1) then
	 a=a(2:)
	 goto 10
      else
	 b=a(:i-1)
	 a=a(i+1:)
      endif
      end

**STRIP1 -- Strip filename (combination) from string A.
*+
      SUBROUTINE STRIP1 (A, B)
      CHARACTER*(*) A, B
*
* This routine takes string A and strips of the first word, or combination
* of words when they are separated by ' - ' or ' + '.
* This is then stored in B while the remaining string is stored back in B.
* If A is empty, A and B will be returned unaltered.
* If A contains one string, B will be A and A will be empty.
*
* Arguments:
*  A  (input): String of words
*    (output): String of words minus the first word(s)
*  B (output): First word(s) in A
*
* Example:
* A = 'The big brown fox' => STRIP1 => B = 'The', A = 'big brown fox'
* A = 'X - Y = Z' => STRIP1 => B = 'X - Y', A = '= Z'
*-
* Replace ' - ' by '(-)' and ' + ' by '(+)'
*
      integer l1,l2

      l1=index(a,'{')
      l2=index(a,'}')
      if (l1.gt.0 .and. l2.gt.l1) call strip10(a(l1:l2),' ','&')
      call strip10(a,' - ','&-&')
      call strip10(a,' + ','&+&')
      call strip10(a,' * ','&*&')
      call strip(a,b)
      call strip10(a,'&',' ')
      call strip10(b,'&',' ')
      end

      subroutine strip10(a,b,c)
      character*(*) a,b,c
      integer i,l

      l=len(c)-1

10    i=index(a,b)
      if (i.ne.0) then
	 a(i:i+l)=c
	 goto 10
      endif
      end
