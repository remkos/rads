/*MALLOCF/DALLOCF -- Variable memory allocation
*&MALLOCF -- Allocate memory block
*+
      FUNCTION MALLOCF (SIZE, POINTER)
      INTEGER MALLOCF, SIZE, POINTER

  This function allocates a memory block of SIZE bytes. At return,
  POINTER is the pointer to the memory block. If the allocation is
  successful, MALLOCF will be a zero.

  The function should be called as:

      IER = MALLOCF (SIZE, POINTER)

  Then the memory block can be used in your program. Obviously, you can
  not use POINTER as an array. You should use it in this fashion, e.g.,

      CALL FILL (SIZE/4, %val(POINTER))

  with the routine FILL defined as:

      SUBROUTINE FILL(N,X)
      INTEGER N,I
      REAL X(*)
      DO I=1,N
         X(I)=1e0
      ENDDO
      END

  To deallocate the memory (i.e., to give the memory free), use
  the function DALLOCF

  Arguments:
   SIZE     (input): Size of the memory block in bytes
   POINTER (output): Pointer to the memory block
   MALLOCF (output): = 0 on proper return, = 1 on error.
--
*&DALLOCF -- Deallocate memory block
*+
      SUBROUTINE DALLOCF (POINTER)
      INTEGER POINTER

  This subroutine deallocates a memory block. This means
  that the memory allocated by MALLOCF will be given free for other
  use. POINTER is the pointer to the memory block.

  The routine should be called as:

       CALL DALLOCF (POINTER)

  Then the memory block can be used again by others.

  Argument:
   POINTER  (input): Pointer to the memory block
--
  22-Jul-2005 - Made part of RADS
   1-Aug-2005 - Turned DALLOCF into routine and remove SIZE argument
------------------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>

#include <sysdep.h>
#ifdef CAPITALS
#define mallocf MALLOCF
#define dallocf DALLOCF
#endif
#ifdef UNDERSCOREBEFORE
#define mallocf _mallocf
#define dallocf _dallocf
#endif
#ifdef UNDERSCOREAFTER
#define mallocf mallocf_
#define dallocf dallocf_
#endif

int mallocf(size, pointer)
int *size;
void **pointer;
{
    *pointer = malloc(*size);
    if (*pointer == NULL) return 1;
    return 0;
}

int dallocf(pointer)
void **pointer;
{
    free(*pointer);
    return 0;
}
