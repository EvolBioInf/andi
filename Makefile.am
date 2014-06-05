SUBDIRS = . libs/RMQ src docs

noinst_DATA = README

README: README.md
	cp README.md README
	
README.md:
	cp README README.md

.PHONY: code-docs
code-docs:
	doxygen

