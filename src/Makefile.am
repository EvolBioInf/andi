bin_PROGRAMS = andi
andi_SOURCES = andi.cpp esa.cpp process.cpp sequence.cpp io.cpp global.h esa.h process.h sequence.h io.h dist_hack.h $(top_srcdir)/libs/kseq.h
andi_CPPFLAGS = $(OPENMP_CXXFLAGS) -std=gnu++0x -fpermissive -I$(top_srcdir)/libs -I$(top_srcdir)/libs/RMQ -I$(top_srcdir)/libs/libdivsufsort-2.0.1/include 
andi_CXXFLAGS = $(OPENMP_CXXFLAGS) -W -Wall
andi_LDADD = $(top_srcdir)/libs/RMQ/libRMQ_succinct.a -L$(top_srcdir)/libs/libdivsufsort-2.0.1/lib/libdivsufsort.la -ldivsufsort

.PHONY: perf
perf: CXXFLAGS= -g -pg -O3 -ggdb
perf: andi
