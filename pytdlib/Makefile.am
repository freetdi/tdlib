include $(top_srcdir)/cython.am

instdir = $(pythondir)/tdlib
inst_PYTHON = __init__.py
inst_LTLIBRARIES = tdlib.la
tdlib_la_SOURCES = tdlib.pyx python_tdlib.cpp

PACKAGE_STRING+= python interface

CLEANFILES += tdlib.cpp

AM_CXXFLAGS = -pthread -Wl,-O1 -Wl,-Bsymbolic-functions \
              -fno-strict-aliasing -g -fwrapv -O2 -Wall \
              -Wformat -Werror=format-security

AM_CPPFLAGS = @PYTHON_INCLUDES@ -I$(top_srcdir)/src
AM_LDFLAGS = -module -omit-version

TESTS = \
	test_approx.py \
	test_apps.py \
	test_bmd.py \
	test_exact.py \
	test_fill_in.py \
	test_lb.py \
	test_misc.py \
	test_post.py \
	test_pre.py \
	test_sep.py

$(TESTS:%.py=%.log): all-am

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_TESTS_ENVIRONMENT = \
    export PYTHONPATH='.';

EXTRA_DIST = \
	$(TESTS) \
	CFGs.py \
	Decomp.py \
	dimacs.py \
	graphs.py \
	Graph.py \
	python_tdlib.hpp \
	RandomGNM_250_1000.py \
	tdlib.pxd \
	tdlib.pyx \
	util.py

$(TESTS:%.py=%.log): P
.PHONY: P
