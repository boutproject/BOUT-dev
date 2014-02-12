# Makefile for the reference and user manuals

.PHONY:all

all: elm_pb.pdf

%.pdf: %.tex
	pdflatex $(@F:.pdf=)
	pdflatex $(@F:.pdf=)

.PHONY:clean

clean: 
	@rm -f *.aux *.log *.bbl *.blg *.toc
	@rm -f elm_reduced.pdf

