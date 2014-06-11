bin_PROGRAMS = andi
andi_SOURCES = andi.cpp esa.cpp process.cpp sequence.cpp global.h esa.h process.h sequence.h $(top_srcdir)/libs/kseq.h
andi_CPPFLAGS = -fopenmp -std=gnu++0x -fpermissive  -I$(top_srcdir)/libs -I$(top_srcdir)/libs/libdivsufsort-2.0.1/include -I$(top_srcdir)/libs/RMQ
andi_CXXFLAGS = -fopenmp
andi_LDADD = $(top_srcdir)/libs/RMQ/libRMQ_succinct.a -L$(top_srcdir)/libs/libdivsufsort-2.0.1/lib/libdivsufsort.la -ldivsufsort

