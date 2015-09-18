ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4
AM_DISTCHECK_CONFIGURE_FLAGS="--enable-unit-tests"

.PHONY: all

if !BUILD_WITH_LIBDIVSUFSORT
OPT_PSUFSORT = opt/psufsort
endif

SUBDIRS = . $(OPT_PSUFSORT) libs opt src docs
DIST_SUBDIRS = . libs opt src docs test opt/psufsort

# Conditionaly build the tests
if BUILD_TESTS

SUBDIRS+= test

TESTS = test/test_esa test/test_seq test/test_random.sh test/test_join.sh

$(TESTS): src/andi
#	$(MAKE) -C test

endif # BUILD_TESTS


dist_noinst_DATA = README ChangeLog README.md andi-manual.pdf

# This is a small hack to satisfy both GitHub and the GNU coding style.
README: README.md
	cp README.md README


# Recreate the changelog, when the version string changes.
ChangeLog: configure.ac
	echo "Missing Git" > ChangeLog;
	if test -d $(srcdir)/.git; then \
		which git && git log --stat --date=short --abbrev-commit | grep --invert-match '^ [[:alnum:].]' | git stripspace > ChangeLog; \
	fi

.PHONY: code-docs
code-docs:
	cd docs && $(MAKE) code-docs;

docs/manual/andi-manual.pdf: $(SUBDIRS)

# copy the manual to top dir, so it is more visible
andi-manual.pdf: docs/manual/andi-manual.pdf
	cp $^ $@
