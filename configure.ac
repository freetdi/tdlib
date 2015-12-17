#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([tdlib], [0.5.4], [larisch@informatik.uni-frankfurt.de])
AC_CONFIG_SRCDIR([src/TD_misc.hpp])
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([])
LT_INIT()

# Checks for programs.
AC_PROG_CXX

# python stuff
AC_ARG_WITH([python],
 [AS_HELP_STRING([--with-python],
	[Enable python library])],
 [python=$withval],
 [python='no'])

AS_IF([test "$python" != no], [
python_version_required="${python##python}"
AM_PATH_PYTHON([$python_version_required],,
   AC_MSG_ERROR([python $PYTHON_VERSION is not $python_version_required])
	)

AC_PATH_PROGS([CYTHON],
              [cython2.7 cython],
              [no])
AS_IF([test "$CYTHON" = no], [AC_MSG_ERROR([cannot find cython.])])

PYTHON_INCLUDE_DIR=$PYTHON_PREFIX/include
PYTHON_LIB_DIR=$PYTHON_PREFIX/lib

AC_ARG_WITH(python-include-dir,
            [  --with-python-include-dir=DIR look in DIR for Python headers],
            [PYTHON_INCLUDE_DIR=$withval],)
AC_SUBST(PYTHON_INCLUDE_DIR)

AC_ARG_WITH(python-lib-dir,
            [  --with-python-lib-dir=DIR look in DIR for Python libraries],
            [PYTHON_LIB_DIR=$withval],)
AC_SUBST(PYTHON_LIB_DIR)
])

AM_CONDITIONAL([USE_PYTHON], [test x$python != xno] )

# Checks for libraries.

# Checks for header files.
AS_IF([test x"$PYTHON" != x], [
  AM_CHECK_PYTHON_HEADERS([:],
	[AC_MSG_ERROR(please install python headers matching $PYTHON)])
])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([pytdlib/Makefile])
AC_OUTPUT
