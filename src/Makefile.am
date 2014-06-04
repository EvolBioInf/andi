bin_PROGRAMS = andi
andi_SOURCES = andi.cpp esa.cpp process.cpp sequence.cpp
andi_CPPFLAGS = -fopenmp -std=gnu++0x -fpermissive  -I$(top_srcdir)/libs -I$(top_srcdir)/libs/RMQ
andi_CXXFLAGS = -fopenmp
andi_LDADD = $(top_srcdir)/libs/RMQ/libRMQ_succinct.a

