CC=clang
CPP=g++
VERSION= \"0.3\"

OBJECTS= np.o esa.o process.o sequence.o

INCLUDES= -I ../libs -I../libs/RMQ

FLAGS= -O3 -W -Wall -fpermissive -fopenmp -std=gnu++0x

LIBS= -ldivsufsort -ldivsufsort64 -lgsl -lgslcblas 

PREPROCESS= -D VERSION=$(VERSION)

np: $(OBJECTS)
	$(CPP) $(FLAGS) -o np $(OBJECTS) ../libs/RMQ/RMQ_1_n.o ../libs/RMQ/RMQ_succinct.o -L/usr/local/lib $(LIBS)

%.o: %.c
	$(CPP) $(FLAGS) $(INCLUDES) $(PREPROCESS) $< -c -o $@
	
.PHONY: profile
profile: FLAGS= -W -Wall -fpermissive -std=gnu++0x -g -pg -O1
profile: np

.PHONY : clean
clean:
	-rm *.o

