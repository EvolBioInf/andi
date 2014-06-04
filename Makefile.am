SUBDIRS = . libs/RMQ src

src_DATA = README

README: README.md
	cp README.md README

.PHONY: code-docs
code-docs:
	doxygen

