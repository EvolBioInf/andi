CC=clang
CPP=g++
VERSION= \"pre0.4\"

OBJECTS= andi.o esa.o process.o sequence.o

INCLUDES= -I ../libs -I../libs/RMQ

FLAGS= -O3 -W -Wall -fpermissive -fopenmp -std=gnu++0x

LIBS= -ldivsufsort -lgsl -lgslcblas 

PREPROCESS= -D VERSION=$(VERSION)

andi: $(OBJECTS)
	$(CPP) $(FLAGS) -o andi $(OBJECTS) ../libs/RMQ/RMQ_1_n.o ../libs/RMQ/RMQ_succinct.o -L/usr/local/lib $(LIBS)

%.o: %.c
	$(CPP) $(FLAGS) $(INCLUDES) $(PREPROCESS) $< -c -o $@
	
.PHONY: profile
profile: FLAGS= -W -Wall -fpermissive -std=gnu++0x -g -pg -O1
profile: andi

.PHONY : clean
clean:
	-rm *.o

