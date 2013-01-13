#-----------------------------------------------------------------------
# $Id$
#
# Copyright (c) 2011-2013  Remko Scharroo (Altimetrics LLC)
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

include config.mk

CONFIG_IN   = config.mk.in src/rads-config.in
CONFIG      = $(CONFIG_IN:.in=)
MAKE_IN_DIRS= src $(DEVEL)

#-----------------------------------------------------------------------
# Rules for making and installing
#-----------------------------------------------------------------------

all install clean spotless::	$(CONFIG)
	@for dir in $(MAKE_IN_DIRS); do (cd src; $(MAKE) $@); done

install::
	$(INSTALL_DIR) $(datadir)/conf
	$(INSTALL_DATA) conf/*.xml $(datadir)/conf

clean spotless::
	$(RM) $(CONFIG)

#-----------------------------------------------------------------------
# To make configuration files
#-----------------------------------------------------------------------

config:	$(CONFIG)
$(CONFIG):	configure $(CONFIG_IN)
	@if test -f config.log ; then \
		cmd=`grep '^  \\$$' config.log | grep configure | cut -c5-` ; \
		echo "Running ... $$cmd" ; $$cmd ; \
	else \
		echo "Run configure first" ; \
	fi

configure:	configure.ac
	autoconf

#-----------------------------------------------------------------------
# To make tarballs
#-----------------------------------------------------------------------

tar:
