.PHONY:all
all: gridgen.pdf

%.pdf: %.tex
	latexmk -pdf $< -interaction=batchmode

clean:
	latexmk -C
