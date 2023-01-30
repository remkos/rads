
#ifndef ccf_festide_h
#define ccf_festide_h

void fes_tide_
	(char	*model_name,
         char	*wave_path_name,
         float	*tide,
         float	*tide_lp,
         float	*tideload,
         float	*lat,
         float	*lon,
         long	*jnasa,
         double	*heure,
         long	*iascii,
         long	*istat,
 /* Rajoute automatiquement par le fortran */
 	 long	model_name_length,
         long	wave_path_name_length);

#endif
