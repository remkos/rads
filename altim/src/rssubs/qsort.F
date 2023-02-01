#ifdef DOUBLE
#define QSORT DQSORT
#define TYPE REAL*8
#else
#ifdef SINGLE
#define QSORT RQSORT
#define TYPE REAL*4
#else
#define QSORT IQSORT
#define TYPE INTEGER*4
#endif
#endif
**QSORT -- Make a quick sort of an TYPE array
*+
      SUBROUTINE QSORT (INDX, ARR, NR)
      INTEGER*4 NR, INDX(NR)
      TYPE ARR(*)
*
* This routine quick-sorts an INTEGER*4 array INDX according to the
* corresponding value in array ARR of type TYPE.
* The array INDX contains NR pointers to the values in array ARR, which
* does not have to be equally large.
*
* Arguments:
*  INDX (input) : One-dimensional array of indx numbers (INTEGER*4)
*      (output) : Row of indx numbers sorted on corresponding value
*  ARR  (input) : One-dimensional array of values (TYPE)
*  NR   (input) : Number of indices/values
*
* Example 1:
* Assume you have 6 values 100, 10, 11, 21, 17, and 90, stored in array
* ARR. These values have to be sorted in ascending order. Your input
* will be like this:
*
*    INDX    1    2    3    4    5    6
*     ARR  100   10   11   21   17   90
*      NR    6
*
* After running QSORT, the output will look like this:
*
*    INDX    2    3    5    4    6    1
*     ARR  100   10   11   21   17   90
*      NR    6
*
* Note that the array ARR has not been changed, but that INDX is sorted by
* the order in which ARR should be sorted. You may print the values in
* the ascending order as follows:
*
*     DO I=1,NR
*        WRITE (*,*) ARR(INDX(I))
*     ENDDO
*
* Example 2:
* It is also possible that the indices are not in logical order. As long as
* they point to the respective locations of the values in array ARR, the
* sorting will be done accordingly. Assume you have the following input:
* (values ** are irrelevant)
*
*    INDX    1    2    4    5    6    8
*     ARR  100   10   **   21   17   90   **   1
*      NR    6
*
* Then after running QSORT, the result will be:
*
*    INDX    8    2    5    4    6    1
*     ARR  100   10   **   21   17   90   **   1
*      NR    6
*
* Printing the values in ascending order after this goes the same as in
* Example 1.
*
* This routine is based on the routine indexx() discussed in
* Numerical Recipes in Fortran, The Art of Scientific Computing.
*-
* 30-Mar-1994 - Created by Marc Naeije
* 23-Jun-1994 - Initialisation of INDX removed [RS]
*  4-Aug-1998 - New manual. RQSORT added.
* 10-Aug-1998 - NSTACK increased from 50 to 100
* 16-Jul-1999 - NSTACK increased from 100 to 200
* 12-Jan-2000 - NSTACK increased from 200 to 800
* 27-Aug-2002 - NSTACK increased from 800 to 2000 [ED]
* 22-Jan-2014 - Several bug fixes; NSTACK reduced back to 256 [RS]
*-----------------------------------------------------------------------
      INTEGER*4 M,NSTACK
      PARAMETER (M=7,NSTACK=256)
      INTEGER*4 i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      TYPE a
*     do j=1,nr
*        indx(j)=j
*     enddo
      jstack=0
      l=1
      ir=nr
    1 if (ir-l.lt.M) then
         do j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do i=j-1,l,-1
               if (arr(indx(i)).le.a) goto 2
               indx(i+1)=indx(i)
            enddo
            i=l-1
    2       indx(i+1)=indxt
         enddo
         if (jstack.eq.0) return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         i=l+1
         j=ir
         itemp=indx(k)
         indx(k)=indx(i)
         indx(i)=itemp
         if (arr(indx(l)).gt.arr(indx(ir))) then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if (arr(indx(i)).gt.arr(indx(ir))) then
             itemp=indx(i)
             indx(i)=indx(ir)
             indx(ir)=itemp
         endif
         if (arr(indx(l)).gt.arr(indx(i))) then
            itemp=indx(l)
            indx(l)=indx(i)
            indx(i)=itemp
         endif
         indxt=indx(i)
         a=arr(indxt)
    3    continue
            i=i+1
         if (arr(indx(i)).lt.a) goto 3
    4    continue
            j=j-1
         if (arr(indx(j)).gt.a) goto 4
         if (j.lt.i) goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
    5    indx(l+1)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if (jstack.gt.NSTACK) stop 'NSTACK too small in QSORT'
         if (ir-i+1.ge.j-1) then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      END
