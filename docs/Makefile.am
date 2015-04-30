dist_man_MANS = andi.1
dist_pdf_DATA = Doxyfile manual/andi-manual.pdf

.PHONY: code-docs
code-docs:
	doxygen

manual_files=manual/andi-manual.tex manual/andi_labels.pdf 

manual/andi-manual.pdf: $(manual_files)
	cd manual && pdflatex andi-manual.tex

dist_noinst_DATA= $(manual_files)

maintainer-clean-local:
	rm -f manual/*{aux,log,out,toc} manual/andi-manual.pdf
