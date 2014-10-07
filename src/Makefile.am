bin_PROGRAMS = andi
andi_SOURCES = andi.cpp esa.cpp process.cpp sequence.cpp io.cpp global.h esa.h process.h sequence.h io.h $(top_srcdir)/libs/kseq.h
andi_CPPFLAGS = $(OPENMP_CXXFLAGS) -std=gnu++0x -fpermissive -I$(top_srcdir)/libs
andi_CXXFLAGS = $(OPENMP_CXXFLAGS) -W -Wall

.PHONY: perf
perf: CXXFLAGS= -g -pg -O3 -ggdb
perf: andi
