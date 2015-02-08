AC_INIT([andi], [0.8.2-beta])
AM_INIT_AUTOMAKE([-Wall foreign ])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_LANG(C++)

AC_CHECK_LIB(m, cos)
AC_OPENMP

# For suffix array computation we use libdivsufsort.
# Fail if the library or its headers were not found.
have_libdivsufsort=yes
AC_CHECK_LIB(divsufsort, divsufsort, [], [have_libdivsufsort=no])

if test "x${have_libdivsufsort}" = xyes; then
	AC_CHECK_HEADERS([divsufsort.h],[],[have_libdivsufsort=no])
fi

if test "x${have_libdivsufsort}" = xno; then
	AC_MSG_ERROR([
	The libdivsufsort library and header files are required for compilation. Aborting.
	])
fi


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

	AX_CXX_COMPILE_STDCXX_11([],[mandatory])
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
AC_CHECK_FUNCS([strtoul strtod])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 test/Makefile
 libs/RMQ/Makefile
 docs/Makefile
 docs/andi.1
])
AC_OUTPUT

