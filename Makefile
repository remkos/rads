#-----------------------------------------------------------------------
# Copyright (c) 2011-2015  Remko Scharroo
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
#
# This Makefile provides the following targets:
# all      : Compile code
# install  : Install executables and scripts
# clean    : Remove compiled code
# spotless : Remove compiled code and configuration files
# gitclean : Remove all files not in git repository
# test     : Test the code
# update   : Pull update from git repository
# refresh  : Pull update and install if needed
#
#-----------------------------------------------------------------------

include config.mk

CONFIG_IN = config.mk.in src/rads-config.in
CONFIG    = $(CONFIG_IN:.in=)
SUBDIRS   = src conf $(DEVEL)

#-----------------------------------------------------------------------
# Rules for making, installing, and cleaning
#-----------------------------------------------------------------------

all install clean spotless test::	$(CONFIG)
	@for subdir in $(SUBDIRS); do (echo "Making $@ in $$subdir ..." && cd $$subdir && $(MAKE) $@); done

clean spotless::
	$(RM) $(CONFIG)

gitclean:
	git clean -fdx

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

#-----------------------------------------------------------------------
# To update from GitHub
#-----------------------------------------------------------------------

update:
	git pull origin master

refresh:
	git pull origin master | grep -q "Already up-to-date." ; if test $$? -eq 1 ; then $(MAKE) install ; fi
