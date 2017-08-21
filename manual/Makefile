# Makefile for the reference and user manuals
TEXDIR=tex_files
VPATH=$(TEXDIR)
.PHONY:all clean sphinx

all: sphinx

sphinx:
	sphinx-build -b latex sphinx/ build/
	cd build && latexmk -pdf BOUT
	test -e BOUT.pdf || ln -s build/BOUT.pdf .
	@echo "Documentation is available in BOUT.pdf"

old: user_manual.pdf \
	developer_manual.pdf \
	preconditioning.pdf \
	coordinates.pdf \
	BOUT_Gradperp_op.pdf

%.pdf: %.tex references.bib
	@cd $(TEXDIR) &&\
	pdflatex $(@F:.pdf=) &&\
	if [ $(@F:.pdf=) != BOUT_Gradperp_op ]; \
	then bibtex $(@F:.pdf=); \
	fi &&\
	makeindex $(@F:.pdf=) &&\
	pdflatex $(@F:.pdf=) &&\
	pdflatex $(@F:.pdf=) &&\
	mv *.pdf ./..
	@cd $(TEXDIR) &&\
	rm -f *.dvi *.aux *.out *.log *.toc *.idx *.ilg *.ind *.bbl *.blg

clean:
	@echo "Cleaning up..."
	@cd $(TEXDIR) &&\
	rm -f *.pdf *.dvi *.aux *.out *.log *.toc *.idx *.ilg *.ind *.bbl *.blg
	@rm -f *.pdf
