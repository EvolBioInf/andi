SUBDIRS = . libs/RMQ src docs test

TESTS = test/test_esa test/test_seq

$(TESTS): all
#	$(MAKE) -C test

dist_noinst_DATA = README ChangeLog README.md

README: README.md
	cp README.md README
	
ChangeLog: configure.ac
	which git && git log --stat --date=short --abbrev-commit | grep --invert-match '^ [[:alnum:].]' | git stripspace > ChangeLog

.PHONY: code-docs
code-docs:
	cd docs && $(MAKE) code-docs;
