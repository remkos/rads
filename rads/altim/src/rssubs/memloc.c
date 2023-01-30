/*MEMLOC/MEMGET/MEMPUT -- Direct manipulation of memory location
*+

      FUNCTION MEMLOC (VAR)
      INTEGER*4 MEMLOC, VAR

      SUBROUTINE MEMGET (PNT, NR, VAR)
      INTEGER*4 PNT, NR, VAR

      SUBROUTINE MEMPUT (PNT, NR, VAR)
      INTEGER*4 PNT, NR, VAR

* The function MEMLOC returns the location of a variable in memory.
* The returned value in MEMLOC is basically a pointer to
* variable VAR. VAR can be of any type.
*
* The routine MEMGET gets NR bytes out of the memory location pointed
* to by PNT and places them into variable VAR. VAR can be of any type.
*
* The routine MEMPUT copies NR bytes from the variable VAR and place
* them into the memory location pointed to by PNT. VAR can be of any type.
*
* These routine were created to avoid the use of the %VAL
* construct in subroutine calls.
*-
* 26-July-2006 -- Created by Remko Scharroo, Altimetrics LLC
*-----------------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>

#include <sysdep.h>
#ifdef CAPITALS
#define memloc MEMLOC
#define memput MEMPUT
#define memget MEMGET
#endif
#ifdef UNDERSCOREBEFORE
#define memloc _memloc
#define memput _memput
#define memget _memget
#endif
#ifdef UNDERSCOREAFTER
#define memloc memloc_
#define memput memput_
#define memget memget_
#endif

int memloc (int var)
{
    return var;
}

void memput (char **pnt, int *nr, char *var)
{
	int i;
	char *tmp;
	tmp = *pnt;
	for (i = 0; i < *nr; i++) tmp[i] = var[i] ;
}

void memget (char **pnt, int *nr, char *var)
{
	int i;
	char *tmp;
	tmp = *pnt;
	for (i = 0; i < *nr; i++) var[i] = tmp[i] ;
}
