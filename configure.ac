AC_INIT([andi], [0.4])
AM_INIT_AUTOMAKE([-Wall foreign])

AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_RANLIB

AC_LANG(C++)

AC_CHECK_LIB(m, cos)
AC_CHECK_LIB(divsufsort, divsufsort)

AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_SEARCH_LIBS([gsl_sf_lnchoose], gslcblas)
AC_SEARCH_LIBS([gsl_sf_lnchoose], gsl)

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
])
AC_OUTPUT

