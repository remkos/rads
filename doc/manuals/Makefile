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
SUB = sub/rads.tex sub/rads_geo.tex sub/rads_grid.tex sub/rads_math.tex \
	sub/rads_misc.tex sub/rads_netcdf.tex sub/rads_time.tex
TAB = sub/rads_math.tab

#-----------------------------------------------------------------------
# What to create?
#-----------------------------------------------------------------------

all:	$(MAN)

install:	all
	$(INSTALL_DATA) $(MAN) $(docdir)

rads4_user_manual.pdf:	%.pdf:	%.tex $(CONFIG) $(SUB) $(TAB)
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

#-----------------------------------------------------------------------
# How to create LaTeX documentation from source code
#-----------------------------------------------------------------------

$(SUB):	robodoc.rc
sub/%.tex:	../../src/%.f90
	$(ROBODOC) --src ../../src/$*.f90 --doc $(TMPDIR)/$* --singlefile
	sed 's/\\index{unsorted![^}]*}//g' $(TMPDIR)/$*.tex > sub/$*.tex
	$(RM) $(TMPDIR)/$*.tex


sub/%.tab:	../../src/%.f90
	grep "\!\!" ../../src/$*.f90 | tail -n +2 | cut -c3- | sed -e 's/ : / \& /' | awk '{print $$0,"\\\\"}' > $@

#-----------------------------------------------------------------------
# How to cleanup
#-----------------------------------------------------------------------

clean:
	$(RM) *.aux *.blg *.idx *.ilg *.ind *.log *.gz *.toc
spotles:	clean
	$(RM) rads_config.tex
