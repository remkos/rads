      subroutine roundoff(min,max)
      real min,max,i
      i=(max-min)/2.
      i=log10(i)
      i=10.**(nint(i-0.5))
      min=nint(min/i-0.5)*i
      max=nint(max/i+0.5)*i
      end
