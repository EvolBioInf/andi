AC_INIT([andi], [0.15-beta])
AM_INIT_AUTOMAKE([-Wall foreign])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Make sure, also the C++ programs are compiled with OpenMP
AC_LANG(C++)
AC_OPENMP

# Execute all tests using C
AC_LANG(C)
AC_OPENMP

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm], [], [have_gsl=no])
AC_CHECK_LIB([gsl],[gsl_ran_binomial], [], [have_gsl=no])

AS_IF([test "x$have_gsl" = "xno"],[
	AC_MSG_ERROR([Missing the Gnu Scientific Library.])
])

# The libdivsufsort header contains some Microsoft extension making
# compilation fail on certain systems (i.e. OS X). Add the following
# flag so the build runs smoothly.
CPPFLAGS="$CPPFLAGS -fms-extensions"
AC_CHECK_HEADERS([divsufsort.h],[have_libdivsufsort=yes],[have_libdivsufsort=no])
AC_CHECK_LIB(divsufsort, divsufsort, [], [have_libdivsufsort=no])

AS_IF([test "x$have_libdivsufsort" = "xno"],[
	AC_MSG_ERROR([Missing libdivsufsort.])
])


# The unit tests require GLIB2. So by default do not build the test.
# If enabled, check for glib.

AC_ARG_ENABLE([unit-tests],
	[AS_HELP_STRING([--enable-unit-tests],[build unit tests @<:@default: no@:>@])],
	[try_unit_tests=${enableval}],[try_unit_tests=no]
	)

AM_CONDITIONAL([BUILD_TESTS],[test "x${try_unit_tests}" = xyes])

# The user may set a seed for the unit tests, so that builds are reproducible.
# A value of 0 makes the tests random.
AC_ARG_WITH([seed],
	[AS_HELP_STRING([--with-seed=INT],
		[random seed for reproducible builds. @<:@default: 0@:>@])],
	[SEED=$withval],
	[SEED=0])

AC_SUBST([SEED])

AS_IF([test "x${try_unit_tests}" = xyes], [
	have_glib=yes
	PKG_CHECK_MODULES([GLIB], [glib-2.0], [], [have_glib=no])

	if test "x${have_glib}" = xno; then
		AC_MSG_ERROR([Missing Glib 2. Either install it or build without unit tests.])
	fi

	AX_CXX_COMPILE_STDCXX_11([],[mandatory])
])


# Check for various headers including those used by libdivsufsort.
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h stdint.h inttypes.h err.h errno.h fcntl.h])

AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T
AC_HEADER_STDBOOL

# Until someone convinces me otherwise, I will deactivate the macros
# AC_FUNC_MALLOC and AC_FUNC_REALLOC. They only check if `malloc(0)` retuns a
# non-null pointer. This breaks the build on systems using uClibc, including
# my laptop.
# As requesting zero bytes is not useful, and implementation-defined behaviour,
# it should be avoided in the first place. Thus I really don't need these checks.

AC_CHECK_FUNCS([floor pow sqrt strdup strerror])
AC_CHECK_FUNCS([strndup strcasecmp])
AC_CHECK_FUNCS([strchr strrchr strchrnul])
AC_CHECK_FUNCS([strtoul strtod])
AC_CHECK_FUNCS([reallocarray])

AM_CONDITIONAL([HAVE_REALLOCARRAY], [test "x$ac_cv_func_reallocarray" = xyes])
AM_CONDITIONAL([HAVE_STRCHRNUL], [test "x$ac_cv_func_strchrnul" = xyes])

AC_CONFIG_HEADERS([src/config.h:src/config.hin])

AC_CONFIG_FILES([
 Makefile
 docs/andi.1
 docs/Makefile
 libs/Makefile
 opt/Makefile
 src/Makefile
 test/Makefile
])
AC_OUTPUT

