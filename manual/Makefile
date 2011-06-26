# Makefile for the reference and user manuals

.PHONY:all

all: user_manual.pdf reference.pdf developer_manual.pdf preconditioning.pdf coordinates.pdf

%.pdf: %.tex references.bib
	pdflatex $(@F:.pdf=)
	bibtex $(@F:.pdf=)
	makeindex $(@F:.pdf=)
	pdflatex $(@F:.pdf=)

.PHONY:clean

clean: 
	@rm -f *.aux *.log *.bbl *.blg *.toc
	@rm -f reference.pdf user_manual.pdf
