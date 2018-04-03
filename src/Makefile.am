bin_PROGRAMS = andi

andi_SOURCES = andi.c esa.c process.c sequence.c io.c global.h esa.h process.h sequence.h io.h dist_hack.h \
model.h model.c
andi_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/libs -I$(top_srcdir)/opt -std=gnu99
andi_CFLAGS = $(OPENMP_CFLAGS) -Wall -Wextra -Wno-missing-field-initializers
andi_LDADD = $(PSUFSORT) $(top_builddir)/libs/libpfasta.a $(top_builddir)/opt/libcompat.a

.PHONY: perf
perf: CFLAGS+= -g -O3 -ggdb -fno-omit-frame-pointer
perf: andi
