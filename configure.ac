AC_INIT([andi], [0.6])
AM_INIT_AUTOMAKE([-Wall foreign])

AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_RANLIB
AM_PROG_AR

AC_LANG(C++)

AC_CHECK_LIB(m, cos)

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

AC_OPENMP

AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h])

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

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 libs/RMQ/Makefile
 docs/Makefile
 docs/andi.1
])
AC_OUTPUT

