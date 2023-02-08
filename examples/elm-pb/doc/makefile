# Makefile for the reference and user manuals

.PHONY:all

all: elm_pb.pdf

%.pdf: %.tex
	latexmk -pdf $(@F:.pdf=) -interaction=batchmode

.PHONY:clean

clean:
	latexmk -C

