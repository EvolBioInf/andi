bin_PROGRAMS = andi

if !BUILD_WITH_LIBDIVSUFSORT
PSUFSORT=$(top_builddir)/opt/psufsort/libpsufsort.a
endif

andi_SOURCES = andi.cpp esa.cpp process.cpp sequence.cpp io.cpp global.h esa.h process.h sequence.h io.h dist_hack.h $(top_srcdir)/libs/kseq.h
andi_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/libs
andi_CXXFLAGS = $(OPENMP_CXXFLAGS) -W -Wall
andi_LDADD = $(PSUFSORT)

.PHONY: perf
perf: CXXFLAGS= -g -O3 -ggdb
perf: andi
