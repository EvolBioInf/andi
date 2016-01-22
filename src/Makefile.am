bin_PROGRAMS = andi

if !BUILD_WITH_LIBDIVSUFSORT
PSUFSORT=$(top_builddir)/opt/psufsort/libpsufsort.a
# This is a hack to make sure, andi is created with a C++ compiler
DUMMY=dummy.cxx
endif

andi_SOURCES = andi.c esa.c process.c sequence.c io.c global.h esa.h process.h sequence.h io.h dist_hack.h \
model.h model.c
andi_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/libs -I$(top_srcdir)/opt -std=gnu99
andi_CFLAGS = $(OPENMP_CFLAGS) -Wall -Wextra -Wno-missing-field-initializers
andi_CXXFLAGS = $(OPENMP_CXXFLAGS) -Wall -Wextra
andi_LDADD = $(PSUFSORT) $(top_builddir)/libs/libpfasta.a $(top_builddir)/opt/libcompat.a
nodist_EXTRA_andi_SOURCES = $(DUMMY)

.PHONY: perf
perf: CFLAGS+= -g -O3 -ggdb -fno-omit-frame-pointer
perf: andi
