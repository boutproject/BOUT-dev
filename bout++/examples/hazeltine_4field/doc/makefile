
TARGET=hazeltine_4field

$(TARGET).pdf: $(TARGET).tex #$(TARGET).bib
	pdflatex $(TARGET)
	#bibtex $(TARGET)
	#pdflatex $(TARGET)
	pdflatex $(TARGET)

.PHONY: clean
clean:
	rm $(TARGET).pdf $(TARGET).aux $(TARGET).log
