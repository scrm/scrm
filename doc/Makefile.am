#docfiles = manual.pdf
if HAVE_PDFLATEX
doc_DATA = manual.pdf
MANNAME = manual
MANTEXSRC = $(MANNAME).tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf
DOCDIR=$(abs_top_builddir)/doc
CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANNAME).bbl $(MANNAME).blg $(MANAUX)

$(MANPDF): $(srcdir)/$(MANTEXSRC)
	$(PDFLATEX) $(MANTEXSRC)
	$(BIBTEX) $(MANAUX)
	$(PDFLATEX) $(MANTEXSRC)
	$(PDFLATEX) $(MANTEXSRC)
endif
#dist_doc_DATA = $(docfiles)

EXTRA_DIST = manual.pdf 

#SUBDIRS =
