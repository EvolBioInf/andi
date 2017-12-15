check_PROGRAMS = test_esa test_seq test_fasta test_process
dist_noinst_DATA = test_extra.sh test_random.sh test_join.sh nan.sh low_homo.sh

if !BUILD_WITH_LIBDIVSUFSORT
PSUFSORT=$(top_builddir)/opt/psufsort/libpsufsort.a
# This is a hack to make sure, the test_esa executable is created with a C++ compiler
DUMMY=dummy.cxx
endif

test_seq_SOURCES = test_seq.c $(top_srcdir)/src/sequence.c
test_seq_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/opt -DDEBUG -std=gnu99
test_seq_CFLAGS = -Wall -Wextra $(GLIB_CFLAGS) -Wno-missing-field-initializers
test_seq_LDADD = $(GLIB_LIBS) $(top_builddir)/opt/libcompat.a

test_process_SOURCES = test_process.c $(top_srcdir)/src/esa.c $(top_srcdir)/src/io.c $(top_srcdir)/src/model.c $(top_srcdir)/src/process.c $(top_srcdir)/src/sequence.c $(top_srcdir)/src/global.h
test_process_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/src -I$(top_srcdir)/opt -I$(top_srcdir)/libs -DDEBUG -std=gnu99
test_process_CFLAGS = $(OPENMP_CFLAGS) -Wall -Wextra $(GLIB_CFLAGS) -Wno-missing-field-initializers
test_process_LDADD = $(GLIB_LIBS) $(PSUFSORT) $(top_builddir)/opt/libcompat.a $(top_builddir)/libs/libpfasta.a
test_process_CXXFLAGS = $(OPENMP_CXXFLAGS) -Wall -Wextra
nodist_EXTRA_test_process_SOURCES = $(DUMMY)

test_esa_SOURCES = test_esa.c $(top_srcdir)/src/esa.c $(top_srcdir)/src/sequence.c $(top_srcdir)/src/esa.h
test_esa_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/libs -I$(top_srcdir)/opt -I$(top_srcdir)/src -DDEBUG -std=gnu99
test_esa_CFLAGS = $(OPENMP_CFLAGS) -Wall -Wextra $(GLIB_CFLAGS) -Wno-missing-field-initializers
test_esa_LDADD = $(GLIB_LIBS) $(PSUFSORT) $(top_builddir)/opt/libcompat.a
test_esa_CXXFLAGS = $(OPENMP_CXXFLAGS) -Wall -Wextra
nodist_EXTRA_test_esa_SOURCES = $(DUMMY)

test_fasta_SOURCES = test_fasta.cxx

.PHONY: all
all: $(check_PROGRAMS)
