      program info_FES9521

      implicit real*8 (a-h,o-z)
      parameter (N = 12)
      character*20 model_name
      logical in_sol,pesudo,rad
      real*4 tpd,lat,long
      integer jday,jnasa,istat
      dimension rlat(N), rlon(N),time(N),time1(N), u(3,2), v(3,2)

      data rlat /-65.44040d0, -65.44040d0, -65.44040d0, -65.44040d0,
     &           -65.44040d0, -65.44040d0, 23.68425d0, 23.68425d0,
     &            23.68425d0, 23.68425d0, 23.68425d0, 23.68425d0/

      data rlon /277.80185d0, 277.80185d0, 277.80185d0, 60.34591d0,
     |		60.34591d0, 60.34591d0, 277.80185d0, 277.80185d0,
     |		277.80185d0, 60.34591d0,
     &            60.34591d0, 60.34591d0/

      data time /48157.45835d0, 49354.41959d0, 50486.23869d0,
     |		48157.45835d0, 49354.41959d0, 50486.23869d0,
     |		48157.45835d0, 49354.41959d0, 50486.23869d0,
     |		48157.45835d0, 49354.41959d0, 50486.23869d0/


c...Read sample input data file to get time and position. 
c    rlat - north latitude (in degrees, -90 to  90)
c    rlon - east longitude (in degrees,  0  to 360)
c    time - desired time,  Modified Julian Date in days. e.g.,
c     Jan 1, 1992 is 48622.

      open (unit=54,
     |file='/infocd/Models/FES95.2/data/ortho_load_3.0',
     |form='formatted')

      model_name = 'fes95.2'

      do i = 1,N

c...Convert time to seconds after Jan 1, 1985 desired by the model
c      MJD in days is 46066 for Jan 1, 1985 
      time1(i) = (time(i)-46066)*86400.d0 

c...Convert the time to days and hours
      jday = idint(time1(i)/86400.0)
      jnasa = 9862+jday
      dhour = (time1(i)-86400.0*jday)/3600.0

c...Call FES95.2 model to compute the tide height
c       in_sol - logical denoting whether tide data exists
c                at desired location
c       tpd    - pure ocean tide height in cm
c       tld    - long period tide in cm
c       tod    - loading tide in cm
c       tide   - predicted ocean tide in cm

        in_sol = .TRUE.
        lat = real(rlat(i))
        long = real(rlon(i))
        call otide(model_name,tpd,lat,long,jnasa,dhour,istat)

        if (abs(tpd).gt.900.) in_sol = .FALSE.
        if (istat.lt.4)      in_sol = .FALSE.

c...Compute long period ocean tide height
        call lpeqmt(time(i)*86400.d0,rlat(i),tld)

c...Compute loading ocean tide height
        call csrtptide('load_3.0',
     |		rlat(i),rlon(i),time(i)*86400.d0,tod,in_sol,
     &   u,v,pseudo,rad)

      if (in_sol) then
        tide = tpd+tld+tod
        write(6,*)
        write(6,*) rlat(i),rlon(i),time(i),tide,tpd,tld,tod
      else
        write(6,*)
        write(6,*) rlat(i),rlon(i),time(i)
        write(6,*) 'No tide solution avaliable.'
      endif

      enddo

      stop
      end 

c...The output data is:
c
c   -65.44040000000000       277.8018500000000       48157.45835000000    
c    28.26559285554116       28.53289       1.473965861400622    
c   -1.741267140380947    
c
c   -65.44040000000000       277.8018500000000       49354.41959000000    
c    23.12015948376020       27.32919      -2.304463541753345    
c   -1.904564367674928    
c
c   -65.44040000000000       277.8018500000000       50486.23869000000    
c    57.42044983673921       64.18126      -2.805336981603809    
c   -3.955472336930417    
c
c   -65.44040000000000       60.34591000000000       48157.45835000000    
c   -18.78546885355988      -21.24364       1.473965861400622    
c   0.9842042773490725    
c
c   -65.44040000000000       60.34591000000000       49354.41959000000    
c   -15.70566204191514      -13.77583      -2.304463541753345    
c   0.3746298613494407    
c
c   -65.44040000000000       60.34591000000000       50486.23869000000    
c   -33.05150690035970      -31.73273      -2.805336981603809    
c    1.486559039373984    
c
c    23.68425000000000       277.8018500000000       48157.45835000000    
c   -22.19493745603563      -21.56812     -0.5132170289758495    
c  -0.1136032853361518    
c
c    23.68425000000000       277.8018500000000       49354.41959000000    
c   -19.46987177979936      -20.49364      0.8023862446569660    
c   0.2213809678532443    
c
c    23.68425000000000       277.8018500000000       50486.23869000000    
c   -12.03505069969826      -13.69473      0.9767842992012572    
c   0.6828919450702118    
c
c    23.68425000000000       60.34591000000000       48157.45835000000    
c   -42.55495786924797      -43.68769     -0.5132170289758495    
c    1.645953709288426    
c
c    23.68425000000000       60.34591000000000       49354.41959000000    
c    16.01606591848306       15.32460      0.8023862446569660    
c  -0.1109186385518337    
c
c    23.68425000000000       60.34591000000000       50486.23869000000    
c    78.52961138758124       79.95807      0.9767842992012572    
c   -2.405241759276268 


