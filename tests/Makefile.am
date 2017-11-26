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

if USE_GALA
GPROG = \
	gala_bmd \
	ta
else
GPROG =
endif

check_PROGRAMS = \
	${GPROG} \
	bmd \
	directed_view \
	bag \
	bucket_sorter_1 \
	bucket_sorter_2 \
	bucket_sorter_2u \
	common_out_edges \
	custom_td \
	dead_simple \
	elim0 \
	inline \
	graph \
	pp \
	pp0 \
	print_concept \
	sdcc_bb \
	sdcc_ub \
	sdcc_uu \
	sdcc_bu \
	vec2tree

# these fail. known bug...
XFAIL_TESTS = \
	bucket_sorter_2u

inline_SOURCES=sdcc_bb.cpp inline.cpp
inline_CXXFLAGS=${AM_CXXFLAGS}

notyets = \
	svbs_random \
	onion \
	gala_hw \
	dfs \
	bfs

${TESTS:%=%.o}: treedec

TESTS = $(check_PROGRAMS)

AM_DEFAULT_SOURCE_EXT=.cpp

BUILT_SOURCES = treedec
CLEANFILES = ${BUILT_SOURCES} *~

AM_CPPFLAGS = -I. -I${srcdir} -DTRACE_UNTESTED -UNDEBUG
AM_CXXFLAGS = -O0 -g

treedec:
	-ln -sf $(top_srcdir)/src treedec

EXTRA_DIST = g.h sdcc_.cpp
# 	svbs_random.h

${TESTS:%=%.log}: .P

.PHONY: .P

.P:
	@:
