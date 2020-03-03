
TARGET = bout_mdsplus

$(TARGET).pdf: $(TARGET).tex
	latexmk -pdf $(TARGET).tex -interaction=batchmode

clean:
	latexmk -C
