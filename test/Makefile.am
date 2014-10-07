check_PROGRAMS = test_esa test_seq

test_seq_SOURCES = test_seq.cxx ../src/sequence.cpp
test_seq_CPPFLAGS = -std=gnu++0x -I$(top_srcdir)/src -DDEBUG
test_seq_CXXFLAGS = -W -Wall $(GLIB_CFLAGS)
test_seq_LDADD = $(GLIB_LIBS)

test_esa_SOURCES = test_esa.cxx ../src/esa.cpp ../src/sequence.cpp $(top_srcdir)/src/esa.h
test_esa_CPPFLAGS = -std=gnu++0x -I$(top_srcdir)/libs -I$(top_srcdir)/src -DDEBUG
test_esa_CXXFLAGS = -W -Wall $(GLIB_CFLAGS)
test_esa_LDADD = $(GLIB_LIBS)

.PHONY: all
all: CXXFLAGS = $(GLIB_CFLAGS) -g -O0 -W -Wall
all: $(check_PROGRAMS)
