check_PROGRAMS = test_esa test_seq test_fasta
dist_noinst_DATA = test_random.sh test_join.sh

if !BUILD_WITH_LIBDIVSUFSORT
PSUFSORT=$(top_builddir)/opt/psufsort/libpsufsort.a
endif

test_seq_SOURCES = test_seq.c ../src/sequence.c
test_seq_CPPFLAGS = -I$(top_srcdir)/src -DDEBUG -std=gnu11
test_seq_CFLAGS = -W -Wall $(GLIB_CFLAGS)
test_seq_LDADD = $(GLIB_LIBS)

test_esa_SOURCES = test_esa.c ../src/esa.c ../src/sequence.c $(top_srcdir)/src/esa.h
test_esa_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/libs -I$(top_srcdir)/src -DDEBUG -std=gnu11
test_esa_CFLAGS = $(OPENMP_CXXFLAGS) -W -Wall $(GLIB_CFLAGS)
test_esa_LDADD = $(GLIB_LIBS) $(PSUFSORT)

test_fasta_SOURCES = test_fasta.cxx

.PHONY: all
all: $(check_PROGRAMS)
