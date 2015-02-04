AC_INIT([andi], [0.8.1-beta])
AM_INIT_AUTOMAKE([-Wall foreign ])

AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_LANG(C++)

AC_CHECK_LIB(m, cos)
AC_OPENMP

# By default try to build with libdivsufsort.
AC_ARG_WITH([libdivsufsort],
    AS_HELP_STRING([--without-libdivsufsort], [Build without libdivsufsort and use psufsort instead.]))

AS_IF([test "x$with_libdivsufsort" != "xno"],
	[
		AC_CHECK_HEADERS([divsufsort.h],[have_libdivsufsort=yes],[have_libdivsufsort=no])
		AC_CHECK_LIB(divsufsort, divsufsort, [], [have_libdivsufsort=no])
	],
	[have_libdivsufsort=no]
)

AS_IF([test "x$have_libdivsufsort" = "xno"],
	[AS_IF([test "x$with_libdivsufsort" = "xyes"],
		[AC_MSG_ERROR([Configuration for libdivsufsort failed. Either install libdivsufsort,
			or use out replacement psufsort instead (./configure --without-libdivsufsort).
			The latter may result in longer runtimes.])
	])
])

AM_CONDITIONAL([BUILD_WITH_LIBDIVSUFSORT],[test "x${with_libdivsufsort}" = xyes])

# The unit tests require GLIB2. So by default do not build the test.
# If enabled, check for glib.

AC_ARG_ENABLE([unit-tests],
	[AS_HELP_STRING([--enable-unit-tests],[build unit tests @<:@default: no@:>@])],
	[try_unit_tests=${enableval}],[try_unit_tests=no]
	)

AM_CONDITIONAL([BUILD_TESTS],[test "x${try_unit_tests}" = xyes])

if test "x${try_unit_tests}" = xyes; then

	have_glib=yes
	PKG_CHECK_MODULES([GLIB], [glib-2.0], [], [have_glib=no])

	if test "x${have_glib}" = xno; then
		AC_MSG_ERROR([Missing Glib 2. Either install it or build without unit tests.])
	fi

fi

# Check for various headers including those used by libdivsufsort.
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h stdint.h inttypes.h err.h errno.h])

AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T
AC_HEADER_STDBOOL

AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt strdup])
AC_CHECK_FUNCS([strndup])
AC_CHECK_FUNCS([strrchr])

AC_CONFIG_HEADERS([src/config.h:src/config.hin])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 test/Makefile
 libs/RMQ/Makefile
 docs/Makefile
 docs/andi.1
 opt/psufsort/Makefile
])
AC_OUTPUT

