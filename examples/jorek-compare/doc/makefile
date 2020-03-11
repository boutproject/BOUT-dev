
DOCS=jorek_compare.pdf

all:$(DOCS)

%.pdf: %.tex
	latexmk -pdf $< -interaction=batchmode

.PHONY: clean
clean:
	latexmk -C
