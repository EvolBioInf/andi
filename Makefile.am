SUBDIRS = . libs/RMQ src docs

noinst_DATA = README

README: README.md
	cp README.md README

.PHONY: code-docs
code-docs:
	doxygen

