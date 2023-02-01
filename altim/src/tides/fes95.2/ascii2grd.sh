#!/bin/csh -fx
foreach comp ( $* )
  ../ascii2grd $comp.fes95.2
  cmp $comp.fes95.2.amp /home/altim/data/grenoble/$comp.fes95.2.1.amp
  cmp $comp.fes95.2.pha /home/altim/data/grenoble/$comp.fes95.2.1.pha
end
