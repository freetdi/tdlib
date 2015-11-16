include $(top_srcdir)/cython.am

instdir = $(pythondir)/tdlib
inst_PYTHON = __init__.py
inst_LTLIBRARIES = tdlib.la
tdlib_la_SOURCES = tdlib.pyx python_tdlib.cpp

CLEANFILES += tdlib.cpp

AM_CXXFLAGS = -pthread -Wl,-O1 -Wl,-Bsymbolic-functions \
              -fno-strict-aliasing -g -fwrapv -O2 -Wall \
              -Wformat -Werror=format-security

AM_CPPFLAGS = -DNDEBUG @my_CPPFLAGS@ -I$(top_srcdir)/src
AM_LDFLAGS = -module -omit-version

TESTS = test.py
TEST_EXTENSIONS = .py
PY_LOG_COMPILER = python
AM_TESTS_ENVIRONMENT = \
    export PYTHONPATH='.';

EXTRA_DIST = tdlib.pxd tdlib.pyx python_tdlib.hpp $(TESTS)
