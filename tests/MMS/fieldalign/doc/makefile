
DOCS = fieldalign.pdf

.PHONY: all

all: $(DOCS)

%.pdf: %.tex 
	latexmk -pdf $< -interaction=batchmode

clean:
	latexmk -C
