CC=clang
CPP=g++

OBJECTS= np.o esa.o divergence.o process.o sequence.o

INCLUDES= -I ../libs -I../libs/RMQ

FLAGS= -O1 -W -Wall -fpermissive -fopenmp -std=gnu++0x

LIBS= -ldivsufsort64 -lgsl -lgslcblas 

np: $(OBJECTS)
	$(CPP) $(FLAGS) -o np $(OBJECTS) ../libs/RMQ/RMQ_1_n.o ../libs/RMQ/RMQ_succinct.o -L/usr/local/lib $(LIBS)

%.o: %.c
	$(CPP) $(FLAGS) $(INCLUDES) $< -c -o $@
	


.PHONY : clean
clean:
	-rm *.o

