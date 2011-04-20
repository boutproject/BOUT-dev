
DOCS=jorek_compare.pdf

all:$(DOCS)

%.pdf: %.tex 
	pdflatex $(@F:.pdf=.tex)
	pdflatex $(@F:.pdf=.tex)

.PHONY: clean
clean:
	rm $(DOCS)
