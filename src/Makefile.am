bin_PROGRAMS = andi

if !BUILD_WITH_LIBDIVSUFSORT
PSUFSORT=$(top_builddir)/opt/psufsort/libpsufsort.a
# This is a hack to make sure, andi is created with a C++ compiler
DUMMY=dummy.cxx
endif

andi_SOURCES = andi.c esa.c process.c sequence.c io.c global.h esa.h process.h sequence.h io.h dist_hack.h
andi_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/libs -I$(top_srcdir)/opt -std=gnu99
andi_CFLAGS = $(OPENMP_CFLAGS) -W -Wall -Wno-missing-field-initializers
andi_CXXFLAGS = $(OPENMP_CXXFLAGS) -W -Wall
andi_LDADD = $(PSUFSORT) $(top_builddir)/libs/libpfasta.a $(top_builddir)/libcompat.a
nodist_EXTRA_andi_SOURCES = $(DUMMY)

.PHONY: perf
perf: CFLAGS+= -g -O3 -ggdb
perf: andi

.PHONY: asan
asan: CFLAGS+= -fsanitize=address -fno-omit-frame-pointer -g -O1
asan: CXXFLAGS+= -fsanitize=address -fno-omit-frame-pointer -g -O1
asan: andi
