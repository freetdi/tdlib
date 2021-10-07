
# legacy tests, with wrapper...
TDLIB_TESTS = \
	test_long_fi_Dimacs.py \
	test_long_approx_CFGs.py \
	dclc_g5.py \
	exact_ta.py \
	test_pp_fi.py

PY_TESTS = \
	graph_bal0.py \
	graph_gala0.py \
	excut_0.py \
	exta_0.py \
	fi_0.py \
	fi_dimacs.py \
	order0.py \
	pp_dimacs.py \
	pp_cfg.py \
	misc_0.py

TESTS = \
	${TDLIB_TESTS} \
	${PY_TESTS}

# AM_PY_LOG_DRIVER_FLAGS=--short

AM_PY_LOG_DRIVER_FLAGS=$(DEBUGFLAG)

# pretend there is a "tdlib" module
all-am: | treedec

treedec:
	- $(LN_S) $(top_srcdir)/tdlib/graphs.py   #HACK
	- $(LN_S) $(top_srcdir)/tdlib/tests/Dimacs.py   #HACK
	- $(LN_S) $(top_srcdir)/tdlib/tests/CFGs.py   #HACK
	- $(LN_S) $(top_srcdir)/tdlib/tests/NewCFGs.py   #HACK
	$(LN_S) ../ treedec                     # pretend there was a treedec module

$(TESTS:%.py=%.log): all-am

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(top_srcdir)/py_log_compiler
AM_TESTS_ENVIRONMENT = \
    export PYTHON=python3; \
	 export PYTHONPATH='.':${srcdir}:${srcdir}/tdlib;

TDLIB_COMPAT = \
	README \
	tdlib.py \
	base.py \
	Graph.py \
	Decomp.py

EXTRA_DIST = \
	$(PY_TESTS) \
	${TESTS:%=ref/%.out} \
	$(TDLIB_COMPAT:%=tdlib/%)

all-local:
	-${LN_S} ${srcdir}/ref .

clean-local:
	rm -rf __pycache__

SED_OUT = cat

${TDLIB_TESTS:%.py=%.py.out}:

${PY_TESTS:%.py=%.py.out}: %.py.out: %.py
	PYTHONPATH=${PYTHONPATH}:${srcdir}:${srcdir}/tdlib \
	${PYTHON} $< | ${SED_OUT} > $@; exit $${PIPESTATUS[0]}

${TDLIB_TESTS:%.py=%.py.out}: %.py.out: %.py
	PYTHONPATH=${PYTHONPATH}:${srcdir}:${srcdir}/tdlib \
	${PYTHON} $< short | ${SED_OUT} > $@; exit $${PIPESTATUS[0]}


$(TESTS:%.py=%.log): P
$(TESTS:%.py=%.py.out): P
.PHONY: P

VPATH += ${top_srcdir}/tdlib/tests
