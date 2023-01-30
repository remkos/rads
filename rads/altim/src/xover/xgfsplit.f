      program xgfsplit

      character*80 filenm,spec*4
      logical asc1,asc2/.false./,asc(50000)
      integer*4 list(50000),rec(4)
      real*8 t0,rev,dlat,r,x
      integer*4 na,nd,narg,iarg,iargc,i,nrec,lat0,lat1,lat2,nlist,
     |		j,i0,i1

      rev=0
      na=0
      nd=0
      narg=iargc()

      if (narg.lt.3) goto 1300
      call getarg(1,filenm)
      open (10,file=filenm,status='old',form='unformatted',
     .recl=18,access='direct')
      read (10,rec=1) spec,nrec
      if (spec.ne.'@XGF') goto 1310
      call getarg(2,filenm)
      open (21,file=filenm,status='new',form='unformatted',
     .recl=18,access='direct')
      call getarg(3,filenm)
      open (22,file=filenm,status='new',form='unformatted',
     .recl=18,access='direct')
      do iarg=4,narg
	 call getarg(iarg,filenm)
	 if (filenm(1:5).eq.'name=') then
	    read (filenm(6:),*) t0,rev
	 endif
      enddo

      if (rev.ne.0) goto 300
* Compute nodes
      read (10,rec=2) rec
      lat0=rec(2)
      read (10,rec=3) rec
      lat1=rec(2)
      nlist=0
      do i=4,nrec+1
	 read  (10,rec=i) rec
	 lat2=rec(2)
         asc1=(lat1.gt.lat0)
	 asc2=(lat2.gt.lat1)
	 if (asc1.neqv.asc2) then
            call inter(dble(lat0),dble(lat1),dble(lat2),x,dlat)
	    nlist=nlist+1
	    if (x.lt.0) then
	       list(nlist)=i-2
	    else
	       list(nlist)=i-1
	    endif
	    asc(nlist)=asc1
*	    write (6,*) nlist,list(nlist),asc(nlist)
	 endif
	 lat0=lat1
	 lat1=lat2
      enddo
      nlist=nlist+1
      list(nlist)=nrec+1
      asc(nlist)=asc2
*      write (6,*) nlist,list(nlist),asc(nlist)
      i0=2
      do j=1,nlist
	 i1=list(j)
	 do i=i0,i1
	    read (10,rec=i) rec
	    if (asc(j)) then
	       na=na+1
	       write (21,rec=na+1) rec
	    else
	       nd=nd+1
	       write (22,rec=nd+1) rec
	    endif
	 enddo
	 i0=i1+1
      enddo
      write (21,rec=1) spec,na
      write (22,rec=1) spec,nd
      write (6,*) 'total     : ',nrec
      write (6,*) 'ascending : ',na
      write (6,*) 'descending: ',nd
      goto 9999

300   continue
* Determine asc/des from time of node and rev length
      do i=2,nrec+1
	 read (10,rec=i) rec
	 r=(rec(1)-t0)/rev
	 if (mod(nint(r*2),2).eq.0) then
* Ascending
	    na=na+1
	    write (21,rec=na+1) rec
	 else
	    nd=nd+1
	    write (22,rec=nd+1) rec
	 endif
      enddo
      write (21,rec=1) spec,na
      write (22,rec=1) spec,nd
      write (6,*) 'total     : ',nrec
      write (6,*) 'ascending : ',na
      write (6,*) 'descending: ',nd
      goto 9999
            

1300  write (6,1301)
1301  format ('usage: xgfsplit <xgfin> <xgfasc> <xgfdes> [ options ]'//
     .'where ''options'' are:'/
     .7x,'node=t0,rev : Specify time of asc node and length of',
     .' revolution (sec)')
      goto 9999
1310  write (6,*) 'xgfsplit: illegal input file'
9999  end

      subroutine inter(x0,x1,x2,t,xt)
      implicit real*8 (a-h,o-z)
      b=(x2-x0)/2
      c=x2-b-x1
      t=-b/(2*c)
      xt=x1+b*t+c*t*t
      end
