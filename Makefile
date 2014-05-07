.PHONY:all
all: 
	cd ./libs/RMQ; make
	cd ./src; make
	@mv ./src/np .

.PHONY: profile
profile:
	cd ./libs/RMQ; make profile
	cd ./src; make profile
	@mv ./src/np .
	
clean:
	@cd ./libs/RMQ; make clean
	@cd ./src; make clean
	@-rm np
	
