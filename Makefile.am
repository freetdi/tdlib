## Process this file with automake to produce Makefile.in
#
# Makefile.am
#
# Automake input file.
#
# Copyright (C) 2015-2016 Felix Salfelder
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this package; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
# Boston, MA 02110-1301, USA.
#

ACLOCAL_AMFLAGS = -Im4

if USE_PYTHON
if DIST
treedec_dir = treedec
endif
if USE_CYTHON
pytdlib_dir = tdlib
endif
endif

SUBDIRS = tests $(treedec_dir) examples $(pytdlib_dir) src


EXTRA_DIST = \
	log_compiler \
	py_log_compiler \
	cython.am
