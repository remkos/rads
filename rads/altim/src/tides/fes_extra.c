/*
 File :       fes_extra.c
 Brief :      Extras for libfes
 Author :     EUMETSAT
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fes.h"

#include "error.h"
#include "fes_int.h"
#include "ini.h"
#include "grid.h"
#include "prediction.h"

/*
 */
void fes_set_nodal_time(void* handle, double time)
{
    fes_handler* fes = handle;

    if (fes != NULL) fes->nodal_time = time;
}
