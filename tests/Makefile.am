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
	bug_19 \
	bucket_sorter_1 \
	bucket_sorter_2 \
	bucket_sorter_3 \
	common_out_edges \
	custom_td \
	dead_simple \
	elim0 \
	inline \
	iter0 \
	iter1 \
	graph \
	pp \
	pp0 \
	pp1 \
	pp_path \
	print_concept \
	sdcc_bb \
	sdcc_ub \
	sdcc_ub_set \
	sdcc_uu \
	sdcc_bu \
	vec2tree

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
LOG_COMPILER = $(top_srcdir)/log_compiler

AM_DEFAULT_SOURCE_EXT=.cpp

BUILT_SOURCES = treedec
CLEANFILES = ${BUILT_SOURCES} *~

AM_CPPFLAGS = -I. -I${srcdir} -DTRACE_UNTESTED -UNDEBUG
AM_CXXFLAGS = -O0 -g

treedec:
	-ln -sf $(top_srcdir)/src treedec

all-local:
	-${LN_S} ${srcdir}/ref .

EXTRA_DIST = \
	bucket_sorter_2u.cpp \
	g.h sdcc_.cpp g_19.h
# 	svbs_random.h

${TESTS:%=%.log}: .P

.PHONY: .P

.P:
	@:
