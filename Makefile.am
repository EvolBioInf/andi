SUBDIRS = . libs/RMQ libs/libdivsufsort-2.0.1 src docs

noinst_DATA = README

README: README.md
	cp README.md README
	
README.md:
	cp README README.md

.PHONY: code-docs
code-docs:
	doxygen

