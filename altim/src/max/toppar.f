
      function toppar(x0,y0,x1,y1,x2,y2)
      real*8 x0,y0,x1,y1,x2,y2,toppar,d
      d=(y0-y1)*(x1-x2)-(y1-y2)*(x0-x1)
      if (abs(d).lt.1d-10) then
	 toppar=1d10
      else
         toppar=((y0-y1)*(x1*x1-x2*x2)-(y1-y2)*(x0*x0-x1*x1))/d/2
      endif
      end
