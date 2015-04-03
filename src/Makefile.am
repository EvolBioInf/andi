bin_PROGRAMS = andi

if !BUILD_WITH_LIBDIVSUFSORT
PSUFSORT=$(top_builddir)/opt/psufsort/libpsufsort.a
# This is a hack to make sure, andi is created with a C++ compiler
DUMMY=dummy.cxx
endif

andi_SOURCES = andi.c esa.c process.c sequence.c io.c global.h esa.h process.h sequence.h io.h dist_hack.h $(top_srcdir)/libs/kseq.h
andi_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/libs -std=gnu11
andi_CFLAGS = $(OPENMP_CFLAGS) -W -Wall -Wno-missing-field-initializers
andi_CXXFLAGS = $(OPENMP_CXXFLAGS) -W -Wall
andi_LDADD = $(PSUFSORT)
nodist_EXTRA_andi_SOURCES = $(DUMMY)

.PHONY: perf
perf: CFLAGS+= -g -O3 -ggdb
perf: andi
