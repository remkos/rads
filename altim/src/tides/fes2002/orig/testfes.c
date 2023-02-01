#include <stdio.h>
#include "festide.h"

/*########################################################################*/

int main(void)
{

/*------------------------------------------------------------------------
 *
 *   PROGRAM : testfes.c
 * 
 *   DESCRIPTION : test the tide prediction algorithm for FES2002
 *                 C version
 * 
 *   PROGRAMMERS : F. LEFEVRE, F. BRIOL
 *  
 *   DATE : June, 3rd 2003
 *------------------------------------------------------------------------
 *
 * rc      : return code (problem if rc != 0)
 * lat     : latitude
 * lon     : longitude
 * time    : time in CNES Julian days
 * hour    : hour
 * tide    : short tides (semi_diurnal and diurnal tides)
 * load    : loading effects
 * lp      : long period tides
 * 
 * tide+lp         = pure tide (as seen by a tide gauge)
 * tide+lp+load    = geocentric tide ((as seen by a satellite)
 * CNES Julian day = NASA Julian day + 2922
 * CNES Julian day 0 is at midnight between the 31 December 
 *                         and 01 January 1950 AD Gregorian
 *
/*########################################################################*/       

  int		rc	= 0;
  double	tide;
  double	load;
  double	lp;

  double	lon	= -7.688;
  double	lat	= 59.195;
  double	time	= 12053;	/* 1983-01-01 00:00:00.0 */

  int		hour	= 0;
  
  /* Initialize memory for FES algorithms */
  if(initFes("../data", fes2002))
    goto onError;

  printf("%12s %5s %9s %9s %9s %9s %9s %9s %9s\n",
    "JulDay","Hour","Latitude","Longitude",
    "Short_tid","LP_tid","Pure_Tide","Geo_Tide","Rad_Tide");

  for(hour = 0; hour < 24; hour++)
  {
    /* Compute tide */
    if(fesTide(lat, lon, time, &tide, &load, &lp))
      goto onError;

    printf("%12.5f %5d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
      time,
      hour, 
      lon, 
      lat,
      tide,
      lp,
      tide+lp,
      tide+lp+load,
      load);

    time += 1/24.0;
  }

  goto onTerminate;

onError:
  rc = 1;

onTerminate:
  /* Free memory for FES algorithms */
  freeFes();

  return rc;
}
