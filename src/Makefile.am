bin_PROGRAMS = andi

if !BUILD_WITH_LIBDIVSUFSORT
PSUFSORT=$(top_builddir)/opt/psufsort/libpsufsort.a
endif

andi_SOURCES = andi.c esa.c process.c sequence.c io.c global.h esa.h process.h sequence.h io.h dist_hack.h $(top_srcdir)/libs/kseq.h
andi_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/libs  -std=gnu11
andi_CFLAGS = $(OPENMP_CFLAGS) -W -Wall
andi_LDADD = $(PSUFSORT)

.PHONY: perf
perf: CFLAGS= -g -O3 -ggdb
perf: andi
