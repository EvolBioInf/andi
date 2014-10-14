bin_PROGRAMS = andi
andi_SOURCES = andi.cpp esa.cpp process.cpp sequence.cpp io.cpp global.h esa.h process.h sequence.h io.h dist_hack.h $(top_srcdir)/libs/kseq.h
andi_CPPFLAGS = $(OPENMP_CXXFLAGS) -std=gnu++0x -fpermissive -I$(top_srcdir)/libs -I$(top_srcdir)/libs/RMQ
andi_CXXFLAGS = $(OPENMP_CXXFLAGS) -W -Wall
andi_LDADD = $(top_builddir)/libs/RMQ/libRMQ_improved.a

.PHONY: perf
perf: CXXFLAGS= -g -pg -O3 -ggdb
perf: andi
