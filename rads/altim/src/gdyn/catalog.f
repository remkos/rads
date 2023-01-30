**CATALOG -- List pass inventory
*+
      SUBROUTINE CATALOG
*-
* 10-Jan-2000 - Only print stats if n>0
*-----------------------------------------------------------------------
      include 'cstg2gbf.inc'
      integer*4 pnt(maxpass)
      character date0*15,date1*15,type*2
      logical	test
      integer*4 i,j,l,n,sat,sta
      integer*4 pnp,onp,pql,oql,t_pnp/0/,t_onp/0/,t_pql/0/,t_oql/0/
      real*8    t0/0d0/,t1/0d0/,t_t0/1d30/,t_t1/-1d30/,value(maxpass)

      test=.true.
      type='QL'
      do l=1,2

* Initialise pointer arrays for QL/NP passes only

      n=0
      do i=1,npass
         if ((p_ql(i).eqv.test)) then
	    n=n+1
	    pnt(n)=i
	 endif
      enddo

* Sort passes by start time

      call dqsort(pnt,p_t0,n)
      
* Print QL/NP passes

      write (*,605) type,n,type
      do i=1,n
         j=pnt(i)
	 call strf1985(date0,'%y%m%d %H:%M:%S',
     |		nint((p_t0(j)-46066)*86400))
	 call strf1985(date1,'%y%m%d %H:%M:%S',
     |		nint((p_t1(j)-46066)*86400))
         write (*,610) j,date0,date1,p_sta(j),p_sat(j),20,p_nr(j)	 
	 if (p_nr(j).eq.0) write (*,*) i,j,p_t0(j),p_t1(j)
      enddo
      test=.false.
      type='NP'
      enddo

* Initialise pointer array and value array for station-satellite list

      do i=1,npass
         value(i)=p_sta(i)+p_sat(i)/1d7
         pnt(i)=i
      enddo

* Sort passes by station-satellite

      call dqsort(pnt,value,npass)

* Make statistics

      write (*,615)
      n=0
      do i=1,npass
         j=pnt(i)
         if (p_sat(j).ne.sat .or. p_sta(j).ne.sta
     |		  .or. i.eq.1) then
	    if (i.ne.1) then
	       n=n+1
	       call strf1985(date0,'%y%m%d %H:%M:%S',
     |		nint((t0-46066)*86400))
	       call strf1985(date1,'%y%m%d %H:%M:%S',
     |		nint((t1-46066)*86400))
	       write (*,620) n,date0,date1,sta,sat,20,pnp,onp,pql,oql
	    endif
	    onp=0
	    oql=0
	    pnp=0
	    pql=0
	    t0=1d30
	    t1=-1d30
	    sat=p_sat(j)
	    sta=p_sta(j)
	 endif
	 t0=min(t0,p_t0(j))
	 t1=max(t1,p_t1(j))
         t_t0=min(t_t0,p_t0(j))
	 t_t1=max(t_t1,p_t1(j))
	 if (p_ql(j)) then
	    pql=pql+1
	    oql=oql+p_nr(j)
	    t_pql=t_pql+1
	    t_oql=t_oql+p_nr(j)
	 else
	    pnp=pnp+1
	    onp=onp+p_nr(j)
	    t_pnp=t_pnp+1
	    t_onp=t_onp+p_nr(j)
	 endif
      enddo
      n=n+1
      call strf1985(date0,'%y%m%d %H:%M:%S',nint((t0-46066)*86400))
      call strf1985(date1,'%y%m%d %H:%M:%S',nint((t1-46066)*86400))
      if (n.gt.0) write (*,620) n,date0,date1,sta,sat,20,pnp,onp,pql,oql

* Print global statistics

      write (*,625)
      call strf1985(date0,'%y%m%d %H:%M:%S',nint((t_t0-46066)*86400))
      call strf1985(date1,'%y%m%d %H:%M:%S',nint((t_t1-46066)*86400))
      if (n.gt.0) write (*,630) n,date0,date1,20,t_pnp,t_onp,t_pql,t_oql
      	 
605   format (/'*** ',a2,' Pass Statistics ***',i7,1x,a2,' Passes ***'/
     |'Passnr            Start               End  Stat',
     |'     SatID Type #Obs')
610   format (i6,2x,a15,' - ',a15,i6,3x,i7.7,i5,i4,2x,a4)
615   format (/'*** Station/satellite Statistics ***'/
     |'     #            First              Last  Stat',
     |'     SatiD Type  #Pas(NP)#Obs  #Pas(QL)#Obs')
620   format (i6,2x,a15,' - ',a15,i6,3x,i7.7,i5,2(i6,i8))
625   format (/'*** Global Statistics ***'/
     |' #Stat            First              Last      ',
     |'           Type  #Pas(NP)#Obs  #Pas(QL)#Obs')
630   format (i6,2x,a15,' - ',a15,6x,3x,7x  ,i5,2(i6,i8))
      end
