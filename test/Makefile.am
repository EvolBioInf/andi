check_PROGRAMS = test_esa
test_esa_SOURCES = test_esa.cxx ../src/esa.cpp ../src/sequence.cpp $(top_srcdir)/src/esa.h
test_esa_CPPFLAGS = -std=gnu++11 -I$(top_srcdir)/libs -I$(top_srcdir)/libs/RMQ -I$(top_srcdir)/src -DDEBUG
test_esa_CXXFLAGS = -W -Wall $(GLIB_CFLAGS)
test_esa_LDADD = $(top_builddir)/libs/RMQ/libRMQ_improved.a $(GLIB_LIBS)

.PHONY: all
all: CXXFLAGS = $(GLIB_CFLAGS) -g -O0 -W -Wall
all: test_esa
