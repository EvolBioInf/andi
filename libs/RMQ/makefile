FLAG = -O3 -fomit-frame-pointer -funroll-loops
CPP=g++

all: testRMQ

testRMQ: testRMQ.cpp RMQ_1_n.o RMQ_succinct.o
	$(CPP) $(FLAG) testRMQ.cpp -o testRMQ  RMQ_1_n.o RMQ_succinct.o

RMQ_succinct.o: RMQ_succinct.cpp RMQ_succinct.hpp RMQ.hpp
	$(CPP) $(FLAG) -c RMQ_succinct.cpp

RMQ_1_n.o: RMQ_1_n.cpp RMQ_1_n.hpp RMQ.hpp
	$(CPP) $(FLAG) -c RMQ_1_n.cpp

.PHONY: profile
profile: FLAG = -O1 -g -pg -fomit-frame-pointer -funroll-loops
profile: all
	

