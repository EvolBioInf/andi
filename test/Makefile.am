check_PROGRAMS = test_esa
test_esa_SOURCES = test_esa.cxx $(top_srcdir)/src/esa.cpp $(top_srcdir)/src/sequence.cpp
test_esa_CPPFLAGS = -std=gnu++0x -I$(top_srcdir)/libs -I$(top_srcdir)/libs/RMQ -DUNIT_TEST -I$(top_srcdir)/src -DDEBUG
test_esa_CXXFLAGS = -W -Wall -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include
test_esa_LDADD = $(top_builddir)/libs/RMQ/libRMQ_improved.a -lglib-2.0

.PHONY: all
all: CXXFLAGS = -g -O0 -W -Wall -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include
all: test_esa
