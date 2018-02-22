
TARGET=drift_instability

$(TARGET).pdf: $(TARGET).tex $(TARGET).bib
	pdflatex $(TARGET)
	bibtex $(TARGET)
	pdflatex $(TARGET)
	pdflatex $(TARGET)

.PHONY: clean
clean:
	rm $(TARGET).pdf
