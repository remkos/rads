        implicit none
        character*20 model
        real*4	 tide, long, lat
        real*8	heure
        integer istat, jnasa
        model='fes95.2.1'
        lat=30.0678
        long=-42.0123 
        jnasa=14230
        heure=4.2343
        do jnasa=14230,14250
        call otide(model,tide,lat,long,jnasa,heure,istat)
        print *,jnasa,tide
        enddo
        stop
        end

