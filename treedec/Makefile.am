AM_DEFAULT_SOURCE_EXT = .cpp
instdir = $(pythondir)/treedec

# TODO: ask python (what to ask?!)
BOOST_PYTHON_LIB_SUFFIX = 311

include MakeList

SUBDIRS = . tests

BOOST_PYTHON_LIBS = \
	-lboost_python${BOOST_PYTHON_LIB_SUFFIX}
	-lboost_numpy${BOOST_PYTHON_LIB_SUFFIX}

_comb_la_LIBADD = ${BOOST_PYTHON_LIBS}
_exact_la_LIBADD = ${BOOST_PYTHON_LIBS}
_graph_la_LIBADD = ${BOOST_PYTHON_LIBS}
_greedy_la_LIBADD = ${BOOST_PYTHON_LIBS}
_lb_la_LIBADD = ${BOOST_PYTHON_LIBS}
_misc_la_LIBADD = ${BOOST_PYTHON_LIBS}
_order_la_LIBADD = ${BOOST_PYTHON_LIBS}
_pp_la_LIBADD = ${BOOST_PYTHON_LIBS}
_treedec_la_LIBADD = ${BOOST_PYTHON_LIBS}

# BUG, not here.
AM_CXXFLAGS = -pthread -Wl,-O1 -Wl,-Bsymbolic-functions \
              -fno-strict-aliasing -g -fwrapv -O2 -Wall \
              -Wformat -Werror=format-security

AM_CPPFLAGS = @PYTHON_INCLUDES@ -I$(top_srcdir)/src
AM_LDFLAGS = -module -avoid-version

EXTRA_DIST = \
	algo_py.hpp \
	boost_compat.h \
	boost_graph.hpp \
	gala_graph_py.hpp \
	graph_iter.hpp \
	graph_py.hpp \
	util_py.hpp

clean-local:
	rm -rf __pycache__


.PHONY: clean-extra symlinks

if VPATH_BUILD
$(PWD)/%.py: | %.py
	$(LN_S) $|
endif

all-local: ${inst_PYTHON:%=$(PWD)/%}
	for i in ${inst_LTLIBRARIES:%.la=%.so}; do \
		rm -f $$i; \
		${LN_S} .libs/$$i $$i; \
	done
