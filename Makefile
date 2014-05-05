.PHONY:all
all: 
	cd ./libs/RMQ; make
	cd ./src; make
	@mv ./src/np .

clean:
	@cd ./src; make clean
	@-rm np
