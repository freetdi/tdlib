#     -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([tdlib], [0.8.4], [larisch@informatik.uni-frankfurt.de])
AC_CONFIG_SRCDIR([src/TD_misc.hpp])
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([])
LT_INIT()

# debugging
AC_MSG_CHECKING([if debug code should be compiled in])
AC_ARG_ENABLE([debug],
[AS_HELP_STRING([--enable-debug],[Enable building of debug code. [[default: disabled]]])],
[
if test "X$enable_debug" = "Xno" ; then
   AC_MSG_RESULT([no])
else
   AC_MSG_RESULT([yes])
   enable_debug=yes
fi
],
[
   AC_MSG_RESULT([no])
   enable_debug=no
])

# Checks for programs.
AC_PROG_CXX

# python stuff
default_python_version=2.7
AC_ARG_WITH([python],
 [AS_HELP_STRING([--with-python],
	[Enable python library])],
 [python=$withval],
 [python="auto"])

AS_IF([test "$python" = auto], [
  AM_PATH_PYTHON([$default_python_version],
					  [python="$default_python_version"], [python="no"]
					 )
],[
AS_IF([test "$python" != no], [
python_version_required="${python##python}"
AM_PATH_PYTHON([$python_version_required],,
   AC_MSG_ERROR([python $PYTHON_VERSION is not $python_version_required])
	)
])])

AS_IF([test "$python" != no], [

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

AS_IF([test "$python" != no], [
	AC_PATH_PROGS([CYTHON], [cython], [no])
	AS_IF([test "$CYTHON" = no], [AC_MSG_ERROR([cannot find cython.])])
])

# Checks for libraries.

# Checks for header files.

# BUG: in auto mode disable python if there are no headers...
AS_IF([test "$python" != no], [
AS_IF([test x"$PYTHON" != x], [
  AM_CHECK_PYTHON_HEADERS([:],
	[AC_MSG_ERROR(please install python headers matching $PYTHON)])
])
])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.

if test "$enable_debug" = "yes" ; then
	CPPFLAGS="$CPPFLAGS -DTRACE_UNTESTED"
else
	CPPFLAGS="$CPPFLAGS -DNDEBUG"
fi

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([pytdlib/Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_OUTPUT

AC_MSG_RESULT([
** Configuration summary for $PACKAGE $VERSION:

   prefix:     $prefix
   python:     $python
   CPPFLAGS:   $CPPFLAGS
])
