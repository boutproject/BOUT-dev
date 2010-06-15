# Makefile for the reference and user manuals

.PHONY:all

all: user_manual.pdf reference.pdf developer_manual.pdf

user_manual.pdf: user_manual.tex references.bib
	pdflatex user_manual
	bibtex user_manual
	pdflatex user_manual
	pdflatex user_manual

reference.pdf: reference.tex
	pdflatex reference.tex
	pdflatex reference.tex

developer_manual.pdf: developer_manual.tex
	pdflatex developer_manual.tex
	pdflatex developer_manual.tex

.PHONY:clean

clean: 
	@rm -f *.aux *.log *.bbl *.blg *.toc
	@rm -f reference.pdf user_manual.pdf
