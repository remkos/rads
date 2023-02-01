      subroutine gridstck(n,m,t,wk1,wk2)
      integer n,m,i,j,k
      real t(m),wk1(n,m),wk2(n,10)
      real*8 sumx,sumy,sumxx,sumxy,sumyy,uxx,uxy,uyy,a,b,r,f,x,y,
     |		miny,maxy

      do i=1,n
         sumx=0d0
         sumy=0d0
         sumxy=0d0
         sumxx=0d0
         sumyy=0d0
	 miny=1d35
	 maxy=-1d35
	 do j=1,m
	    if (wk1(i,j).gt.1e20) then
	       do k=1,8
		  wk2(i,k)=1e30
	       enddo
	       goto 100
	    endif
	    x=t(j)
	    y=wk1(i,j)
            sumx =sumx +x
            sumy =sumy +y
            sumxy=sumxy+x*y
            sumxx=sumxx+x*x
            sumyy=sumyy+y*y
	    miny =min(miny,y)
	    maxy =max(maxy,y)
         enddo
         uxx=m*sumxx-sumx*sumx
         uxy=m*sumxy-sumx*sumy
         uyy=m*sumyy-sumy*sumy
         b=uxy/uxx
	 a=(sumy-b*sumx)/m
         r=uxy/sqrt(uxx*uyy)
	 f=sumyy+b*b*sumxx+m*a*a-2*b*sumxy-2*a*sumy+2*a*b*sumx
	 wk2(i,1)=sumy/m
	 wk2(i,2)=sqrt(sumyy/m)
	 wk2(i,3)=sqrt(sumyy/m-(sumy/m)**2)
	 wk2(i,4)=a
	 wk2(i,5)=b
	 wk2(i,6)=r
	 wk2(i,7)=sqrt(f/m)
	 wk2(i,8)=sumy
	 wk2(i,9)=miny
	 wk2(i,10)=maxy
100      continue
      enddo
      end
