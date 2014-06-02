#OBJECTS= andi.o esa.o process.o sequence.o

#INCLUDES= -I ../libs -I../libs/RMQ

#FLAGS= -O3 -W -Wall -fpermissive -fopenmp -std=gnu++0x

#LIBS= @LIBS@

CC = @CXX@

bin_PROGRAMS = andi
andi_SOURCES = andi.cpp esa.cpp process.cpp sequence.cpp
andi_CPPFLAGS = -fopenmp -std=gnu++0x -fpermissive  -I$(top_srcdir)/libs -I$(top_srcdir)/libs/RMQ
andi_CXXFLAGS = -fopenmp
andi_LDADD = $(top_srcdir)/libs/RMQ/RMQ_succinct.o

#.PHONY: all
#all: andi

#andi: $(OBJECTS)
#	$(CXX) $(FLAGS) -o andi $(OBJECTS) ../libs/RMQ/RMQ_1_n.o ../libs/RMQ/RMQ_succinct.o -L/usr/local/lib $(LIBS)

#%.o: %.c
#	$(CXX) $(FLAGS) $(INCLUDES) $(DEFS) $< -c -o $@
	

