dist_man_MANS = andi.1
dist_pdf_DATA = Doxyfile manual/andi-manual.pdf

.PHONY: code-docs
code-docs:
	doxygen

manual_files=manual/andi-manual.tex manual/img/accuracy/fig-acc2.tex manual/img/accuracy/fig-acc.tex manual/img/concurrency/fig-concurrency.tex manual/img/eco29/resources.tex manual/img/indel/andi.data manual/img/indel/andi_error.data manual/img/indel/cophylog.data manual/img/indel/cophylog_error.data manual/img/indel/fig-error.tex manual/img/indel/fig-indel.tex manual/img/indel/kr.data manual/img/indel/kr_error.data manual/img/indel/spaced.data manual/img/indel/spaced_error.data manual/img/recomb/andi.data manual/img/recomb/cophylog.data manual/img/recomb/fig-recomb.tex manual/img/recomb/kr.data manual/img/recomb/local.data manual/img/recomb/local.tex manual/img/recomb/mugsy.data manual/img/recomb/spaced.data manual/img/threshold/threshold2.data manual/img/threshold/threshold.tex manual/img/3085/andi_sub.pdf manual/img/3085/cophylog_sub.pdf manual/img/concurrency/fig-concurrency.pdf manual/img/eco29/andi_labels.pdf manual/img/eco29/cophylog_labels.pdf manual/img/eco29/kr_labels.pdf manual/img/eco29/mugsy_labels.pdf manual/img/eco29/resources.pdf manual/img/eco29/spaced_k20.pdf manual/img/eco29/spaced_labels.pdf manual/img/indel/fig-error.pdf manual/img/indel/fig-indel.pdf manual/img/recomb/fig-recomb.pdf manual/img/recomb/local.pdf manual/img/threshold/threshold.pdf

manual/andi-manual.pdf: $(manual_files)
	cd manual && pdflatex andi-manual.tex

dist_noinst_DATA= $(manual_files)

maintainer-clean-local:
	rm -f manual/*{aux,log,out,toc} manual/andi-manual.pdf
