ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4
AM_DISTCHECK_CONFIGURE_FLAGS="--enable-unit-tests"

SUBDIRS = . libs/RMQ src docs

# Conditionaly build the tests
if BUILD_TESTS

SUBDIRS+= test

TESTS = test/test_esa test/test_seq test/test_random.sh test/test_join.sh

$(TESTS): all
#	$(MAKE) -C test

endif # BUILD_TESTS


dist_noinst_DATA = README ChangeLog README.md

README: README.md
	cp README.md README
	
ChangeLog: configure.ac
	echo "Missing Git" > ChangeLog;
	if test -d $(srcdir)/.git; then \
		which git && git log --stat --date=short --abbrev-commit | grep --invert-match '^ [[:alnum:].]' | git stripspace > ChangeLog; \
	fi

.PHONY: code-docs
code-docs:
	cd docs && $(MAKE) code-docs;
