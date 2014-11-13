check_PROGRAMS = test_esa test_seq test_fasta

test_seq_SOURCES = test_seq.cxx ../src/sequence.cpp
test_seq_CPPFLAGS = -std=gnu++0x -I$(top_srcdir)/src -DDEBUG
test_seq_CXXFLAGS = -W -Wall $(GLIB_CFLAGS)
test_seq_LDADD = $(GLIB_LIBS)

test_esa_SOURCES = test_esa.cxx ../src/esa.cpp ../src/sequence.cpp $(top_srcdir)/src/esa.h
test_esa_CPPFLAGS = -std=gnu++0x -I$(top_srcdir)/libs -I$(top_srcdir)/libs/RMQ -I$(top_srcdir)/src -DDEBUG -I$(top_srcdir)/libs/libdivsufsort-2.0.1/include 
test_esa_CXXFLAGS = -W -Wall $(GLIB_CFLAGS)
test_esa_LDADD = $(top_builddir)/libs/RMQ/libRMQ_improved.a $(GLIB_LIBS)

test_fasta_SOURCES = test_fasta.cxx
test_fasta_CPPFLAGS = -std=gnu++1y

.PHONY: all
all: CXXFLAGS = $(GLIB_CFLAGS) -g -O0 -W -Wall
all: $(check_PROGRAMS)
