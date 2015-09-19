dist_man_MANS = andi.1
dist_noinst_DATA = Doxyfile

# I intentionally do not list any of the manual files here. I neither want them
# distributed nor installed. The reason is that building the manual requires
# LaTeX with a whole bunch of packages installed. Plus, so many things can go
# wrong, when building, so it's better to inspect the result. Thus, the manual
# has to be build by hand and copied to the right place for distribution.

.PHONY: code-docs
code-docs:
	doxygen

manual/version.tex: manual/version.tex.in $(top_srcdir)/configure.ac
	sed "s/VERSION/$(VERSION)/" manual/version.tex.in > manual/version.tex

manual/andi-manual.pdf: manual/andi-manual.tex manual/version.tex
	@echo "error: manual rebuild of the manual required (no pun intended)."
	@exit 1

# maintainer-clean-local:
#	rm -f manual/*{aux,log,out,toc} manual/andi-manual.pdf
