all: introduction.pdf solvers.pdf code_structure.pdf

%.pdf: %.tex
	latexmk -pdf $< -interaction=batchmode

clean:
	latexmk -C

