
TARGET=drift_instability

$(TARGET).pdf: $(TARGET).tex $(TARGET).bib
	latexmk -pdf $(TARGET) -interaction=batchmode

.PHONY: clean
clean:
	latexmk -C
