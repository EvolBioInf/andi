bin_PROGRAMS = andi
andi_SOURCES = andi.cpp esa.cpp process.cpp sequence.cpp io.cpp global.h esa.h process.h sequence.h io.h dist_hack.h $(top_srcdir)/libs/kseq.h
andi_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/libs
andi_CXXFLAGS = $(OPENMP_CXXFLAGS) -W -Wall

.PHONY: perf
perf: CXXFLAGS= -g -O3 -ggdb
perf: andi
