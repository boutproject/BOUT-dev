
DOCS=dalf3.pdf

all:$(DOCS)

%.pdf: %.tex
	latexmk -pdf $< -interaction=batchmode

.PHONY: clean
clean:
	latexmk -C
