#-----------------------------------------------------------------------
# Copyright (c) 2011-2016  Remko Scharroo
# See LICENSE.TXT file for copying and redistribution conditions.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#-----------------------------------------------------------------------

include ../../config.mk
PDFLATEX = pdflatex
BIBTEX = bibtex
MAKEINDEX = makeindex
ROBODOC = robodoc
CONFIG = rads_config.tex rads4_sats.tex
MAN = rads4_data_manual.pdf rads4_user_manual.pdf
SUBS = rads rads_geo rads_grid rads_math rads_misc rads_netcdf rads_time
SUBS_TEX = $(foreach sub,$(SUBS),sub/$(sub).tex)
PROGS = rads2asc4 rads2nc rads2grd4 radscolin4
PROGS_EXTRA = rads2grd4-o radscolin4-r
PROGS_HLP = $(foreach prog,$(PROGS),sub/$(prog).hlp)
PROGS_OUT = $(foreach prog,$(PROGS) $(PROGS_EXTRA),sub/$(prog).out)
TABS_TEX = sub/rads_math.tab
HLP = sub/rads2asc4.hlp sub
SEPARATOR =	echo "............  lines removed .............................."

#-----------------------------------------------------------------------
# What to create?
#-----------------------------------------------------------------------

all:	$(MAN)

install:	all
	$(INSTALL_DATA) $(MAN) $(docdir)

rads4_user_manual.pdf:	%.pdf:	%.tex $(CONFIG) $(SUBS_TEX) $(PROGS_HLP) $(PROGS_OUT) $(TABS_TEX)
	$(PDFLATEX) $*
	$(PDFLATEX) $*
	$(MAKEINDEX) $*
	$(PDFLATEX) $*

rads4_data_manual.pdf:	%.pdf:	%.tex %.bbl $(CONFIG)
	$(PDFLATEX) $*
	$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*
	$(MAKEINDEX) $*
	$(PDFLATEX) $*

subs:	$(SUBS_TEX) $(PROGS_HLP) $(PROGS_OUT) $(TABS_TEX)

#-----------------------------------------------------------------------
# How to create LaTeX documentation from source code
#-----------------------------------------------------------------------

$(SUBS_TEX):	robodoc.rc
sub/%.tex:	../../src/%.f90
	$(ROBODOC) --src ../../src/$*.f90 --doc $(TMPDIR)/$* --singlefile
	sed 's/\\index{unsorted![^}]*}//g' $(TMPDIR)/$*.tex > sub/$*.tex
	$(RM) $(TMPDIR)/$*.tex

#-----------------------------------------------------------------------
# Make table of the MATH commands
#-----------------------------------------------------------------------

sub/%.tab:	../../src/%.f90
	grep "\!\!" ../../src/$*.f90 | tail -n +2 | cut -c3- | sed -e 's/ : / \& /' | awk '{print $$0,"\\\\"}' > $@

#-----------------------------------------------------------------------
# Make help output
#-----------------------------------------------------------------------

sub/%.hlp:	../../src/%
	PATH=../../src $* --help > $@

#-----------------------------------------------------------------------
# Make command output and logs
#-----------------------------------------------------------------------

sub/rads2asc4.out:	Makefile ../../src/rads2asc4
	PATH=../../src rads2asc4 -S e2 -C 0 -P 901,1000,2 -V time_ymdhms,lat,lon,sla,swh,wind_speed -R -8,42,28,48 -o $(@F) -v > $(@:.out=.log) 0>&2
	(head -n 19 $(@F); $(SEPARATOR); tail -n 4 $(@F)) > $@
	rm -f $(@F)

sub/rads2nc.out:	Makefile ../../src/rads2nc
	PATH=../../src rads2nc -S e2 -C 0 -P 901,1000,2 -V time,lat,lon,sla,swh,wind_speed -R -8,42,28,48 -o $(@F) -v > $(@:.out=.log) 0>&2
	(ncdump $(@F) | head -n 73 | sed 's/	/    /g'; $(SEPARATOR); ncdump $(@F) | tail -n 5) > $@
	rm -f $(@F)

sub/rads2grd4.out:	Makefile ../../src/rads2grd4
	PATH=../../src rads2grd4 -S e2 -C 0 -P 901,1000,2 -V wind_speed,swh,sla -c > $(@F)
	(head -n 16 $(@F); $(SEPARATOR); tail -n 4 $(@F)) > $@
	rm -f $(@F)

sub/rads2grd4-o.out:	Makefile ../../src/rads2grd4
	PATH=../../src rads2grd4 -S e2 -C 0 -P 901,1000,2 -V wind_speed,swh,sla -c -o $(@F)
	(ncdump $(@F) | head -n 42 | sed 's/	/    /g'; $(SEPARATOR); ncdump $(@F) | tail -n 5) > $@
	rm -f $(@F)

sub/radscolin4.out:	Makefile ../../src/radscolin4
	PATH=../../src radscolin4 -S e2 -C 0,2 -P 915,1000,2 -R -8,42,28,48 -V lat,sla,swh > $(@F)
	(head -n 18 $(@F); $(SEPARATOR); tail -n 4 $(@F)) > $@
	rm -f $(@F)

sub/radscolin4-r.out:	Makefile ../../src/radscolin4
	PATH=../../src radscolin4 -S e2 -C 0,2 -P 915,1000,2 -R -8,42,28,48 -V lat,sla,swh -r > $(@F)
	(head -n 18 $(@F); $(SEPARATOR); tail -n 4 $(@F)) > $@
	rm -f $(@F)

inc/max.out:		getraw.nml max.nml $(RADSROOT)/bin/max2
	rm -f xovers.rxf
	max2 sat=e2 cycle=0,1 sel=0,1,17 xovers > inc/max_v.out
	(head -18 xovers.rxf; $(SEPARATOR); tail -4 xovers.rxf) > $@
	rm -f xovers.rxf

#-----------------------------------------------------------------------
# How to cleanup
#-----------------------------------------------------------------------

clean:
	$(RM) *.aux *.blg *.idx *.ilg *.ind *.log *.gz *.toc sub/*
spotless:	clean
	$(RM) rads_config.tex
