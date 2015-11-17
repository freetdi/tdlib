ACLOCAL_AMFLAGS = -Im4

if USE_PYTHON
pytdlib = pytdlib
endif

SUBDIRS = $(pytdlib)

pkginclude_HEADERS = \
  src/TD_simple_graph_algos.hpp \
  src/TD_branch_decomposition.hpp \
  src/TD_combinations.hpp \
  src/TD_dynamicCR.hpp \
  src/TD_elimination_orderings.hpp \
  src/TD_exact_cutset.hpp \
  src/TD_lower_bounds.hpp \
  src/TD_misc.hpp \
  src/TD_NetworkFlow.hpp \
  src/TD_postprocessing.hpp \
  src/TD_preprocessing.hpp \
  src/TD_seperator_algorithm.hpp \
  src/TD_upper_bounds.hpp

EXTRA_DIST = cython.am
