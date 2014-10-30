# Makefile.am for libdivsufsort

noinst_LIBRARIES = libdivsufsort.a
libdivsufsort_a_SOURCES = divsufsort.c sssort.c trsort.c utils.c
libdivsufsort_a_CPPFLAGS = -I$(top_srcdir)/include
#libdivsufsort_a_LDFLAGS = -version-info @LT_CURRENT@:@LT_REVISION@:@LT_AGE@ -export-symbols $(srcdir)/libdivsufsort.sym

if DIVSUFSORT64
noinst_LIBRARIES += libdivsufsort64.a
libdivsufsort64_a_SOURCES = divsufsort.c sssort.c trsort.c utils.c
libdivsufsort64_a_CPPFLAGS = -I$(top_srcdir)/include -DBUILD_DIVSUFSORT64
#libdivsufsort64_a_LDFLAGS = -version-info @LT_CURRENT@:@LT_REVISION@:@LT_AGE@ -export-symbols $(srcdir)/libdivsufsort64.sym
endif

EXTRA_DIST = libdivsufsort.sym libdivsufsort64.sym CMakeLists.txt
