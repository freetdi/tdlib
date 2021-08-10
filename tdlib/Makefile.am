include $(top_srcdir)/cython.am

SUBDIRS = . tests

instdir = $(pythondir)/tdlib
inst_PYTHON = \
	__init__.py \
	graphs.py \
	Decomp.py \
	Graph.py

inst_LTLIBRARIES = cytdlib.la
cytdlib_la_SOURCES = cytdlib.pyx python_tdlib.cpp

PACKAGE_STRING+= python interface

CLEANFILES += cytdlib.cpp

clean-local:
	rm -rf __pycache__

.PHONY: clean-extra symlinks

AM_CXXFLAGS = -pthread -Wl,-O1 -Wl,-Bsymbolic-functions \
              -fno-strict-aliasing -g -fwrapv -O2 -Wall \
              -Wformat -Werror=format-security

longcheck: all-am
	$(MAKE) -C tests $@

if VPATH_BUILD
SYMLINKS = symlinks
endif
all-local: $(SYMLINKS)

symlinks:
	ln -sf $(srcdir)/graphs.py
	ln -sf $(srcdir)/Graph.py
	ln -sf $(srcdir)/Decomp.py
	ln -sf $(srcdir)/printer.py
	ln -sf $(srcdir)/__init__.py

AM_CPPFLAGS = @PYTHON_INCLUDES@ -I$(top_srcdir)/src
AM_LDFLAGS = -module -avoid-version

EXTRA_DIST = \
	graphs.py \
	Graph.py \
	Decomp.py \
	python_tdlib.hpp \
	cytdlib.pxd \
	cytdlib.pyx \
	printer.py

uninstall-hook:
	-rmdir ${DESTDIR}/${instdir}
	-rmdir ${DESTDIR}/${pythondir}
