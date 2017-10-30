
TESTS = $(LONG_TESTS) \
	dclc_g5.py \
	test_generic_search_CFG1.py \
	test_generic_search_CFG2.py \
	test_generic_search_p17.py \
	test_generic_search_gala.py \
	test_approx.py \
	test_apps.py \
	test_bmd.py \
	test_exact.py \
	test_fill_in.py \
	test_lb.py \
	test_misc.py \
	test_post.py \
	test_pp_md.py \
	test_pp_fi.py \
	test_pre.py \
	test_sep.py


LONG_TESTS = \
	test_long_approx_CFGs.py \
	test_long_apps_CFGs.py \
	test_long_lb_CFGs.py \
	test_long_approx_NewCFGs.py \
	test_long_apps_NewCFGs.py \
	test_long_lb_NewCFGs.py \
	test_long_approx_Dimacs.py \
	test_long_lb_Dimacs.py \
	test_long_approx_Zoo.py \
	test_long_apps_Zoo.py \
	test_long_lb_Zoo.py \
        test_long_generic_search_Dimacs.py \
        test_long_generic_search_Networks.py

EXTRALONG_TESTS = \
test_long_fi_Dimacs.py

# AM_PY_LOG_DRIVER_FLAGS=--short

AM_PY_LOG_DRIVER_FLAGS=$(DEBUGFLAG)

longcheck:
	$(MAKE) check TESTS="$(LONG_TESTS)" AM_PY_LOG_DRIVER_FLAGS+=--long

$(TESTS:%.py=%.log): all-am

$(LONG_TESTS:%.py=%.log): AM_PY_LOG_DRIVER_FLAGS=

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_TESTS_ENVIRONMENT = \
    export PYTHONPATH='.:..:$(top_srcdir)/pytdlib';

EXTRA_DIST = \
	$(TESTS) \
	CFGs.py \
	Dimacs.py \
	Zoo.py \
	base.py \
	RandomGNM_250_1000.py \
	util.py

$(TESTS:%.py=%.log): P
.PHONY: P
