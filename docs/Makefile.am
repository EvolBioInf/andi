dist_man_MANS = andi.1
dist_pdf_DATA = manual/andi-manual.pdf

.PHONY: code-docs
code-docs:
	doxygen

manual_files=manual/andi-manual.tex manual/andi_labels.pdf manual/version.tex

manual/version.tex: manual/version.tex.in $(top_srcdir)/configure.ac
	sed "s/VERSION/$(VERSION)/" manual/version.tex.in > manual/version.tex

manual/andi-manual.pdf: $(manual_files)
	-which pdflatex && echo "Missing pdflatex" > $@
	-cd manual && which pdflatex && pdflatex andi-manual.tex

dist_noinst_DATA= Doxyfile $(manual_files)

maintainer-clean-local:
	rm -f manual/*{aux,log,out,toc} manual/andi-manual.pdf
