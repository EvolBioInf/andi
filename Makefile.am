
SUBDIRS = . libs/RMQ src

src_DATA = README

README: README.md
	cp README.md README

#.PHONY:all
#all: 
#	cd ./libs/RMQ; make
#	cd ./src; make
#	@mv ./src/andi .
#
#.PHONY: profile
#profile:
#	cd ./libs/RMQ; make profile
#	cd ./src; make profile
#	@mv ./src/andi .
#
#.PHONY: clean
#clean:
#	@cd ./libs/RMQ; make clean
#	@cd ./src; make clean
#	@-rm andi
#
.PHONY: docs
docs:
	doxygen

