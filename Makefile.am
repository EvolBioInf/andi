SUBDIRS = . libs/RMQ src docs

TESTS = test/test_esa

test/test_esa:
	$(MAKE) -C test

noinst_DATA = README ChangeLog

README: README.md
	cp README.md README
	
README.md:
	cp README README.md

ChangeLog: FORCE
	git log --stat --date=short --abbrev-commit | grep --invert-match '^ [[:alnum:].]' | git stripspace > ChangeLog

.PHONY: code-docs
code-docs:
	cd docs && $(MAKE) code-docs;

.PHONY: FORCE
FORCE:

