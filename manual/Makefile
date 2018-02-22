# Makefile for the reference and user manuals
.PHONY:all clean sphinx html pdf doxygen breathe-autogen

# Some distros may provide sphinx-build-3 and sphinx-build-2 instead
# of sphinx-build. An automatic check would be helpful, but at least
# with this, it should be easy to change:
sphinx-build?=$(shell which sphinx-build-3 &>/dev/null && echo sphinx-build-3 || echo sphinx-build)

all: sphinx
manual: all
# set some shorter names
pdf: sphinx-pdf
html: sphinx-html
man: sphinx-man
sphinx: sphinx-pdf

sphinx-pdf: 
	PYTHONPATH=$(BOUT_TOP)/tools/pylib:$$PYTHONPATH $(sphinx-build) -b latex sphinx/ build/
	cd build && latexmk -pdf BOUT
	test -e BOUT.pdf || ln -s build/BOUT.pdf .
	@echo "Documentation is available in `pwd`/BOUT.pdf"

sphinx-html:
	PYTHONPATH=$(BOUT_TOP)/tools/pylib:$$PYTHONPATH $(sphinx-build) -b html sphinx/ html/
	@echo "Documentation available in `pwd`/html/"

sphinx-man:
	PYTHONPATH=$(BOUT_TOP)/tools/pylib:$$PYTHONPATH $(sphinx-build) -b man sphinx/ man/
	@echo "Documentation available in `pwd`/html/"

# Run doxygen, ignore if it fails (leading '-')
doxygen:
	-cd doxygen && doxygen Doxyfile

# Run breathe-apidoc, ignore if it fails (leading '-')
breathe-autogen: doxygen
	-breathe-apidoc -f -o sphinx/_breathe_autogen doxygen/bout/xml

clean:
	@echo "Cleaning up..."
	@cd $(TEXDIR) &&\
	rm -f *.pdf *.dvi *.aux *.out *.log *.toc *.idx *.ilg *.ind *.bbl *.blg
	@rm -f *.pdf
	@rm -rf html/
	@rm -rf build/
	@rm -rf sphinx/_breathe_autogen/
