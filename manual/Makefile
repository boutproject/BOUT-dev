# Makefile for the reference and user manuals
.PHONY:all clean sphinx html pdf doxygen breathe-autogen

BOUT_TOP?=..

# Some distros may provide sphinx-build-3 and sphinx-build-2 instead
# of sphinx-build. An automatic check would be helpful, but at least
# with this, it should be easy to change:
sphinx-build?=$(shell which sphinx-build-3 >/dev/null 2>&1 && echo sphinx-build-3 || echo sphinx-build)

all: sphinx
manual: all
# set some shorter names
pdf: sphinx-pdf
html: sphinx-html
man: sphinx-man
sphinx: sphinx-pdf

sphinx-pdf: doxygen
	PYTHONPATH=$(BOUT_TOP)/tools/pylib:$$PYTHONPATH $(sphinx-build) -b latex sphinx/ build/
	cd build && latexmk -pdf BOUT -interaction=batchmode
	test -e BOUT.pdf || ln -s build/BOUT.pdf .
	@echo "Documentation is available in $$(pwd)/BOUT.pdf"

sphinx-html: doxygen
	PYTHONPATH=$(BOUT_TOP)/tools/pylib:$$PYTHONPATH $(sphinx-build) -b html sphinx/ html/
	@echo "Documentation available in file://$$(pwd)/html/index.html"

sphinx-man: doxygen
	PYTHONPATH=$(BOUT_TOP)/tools/pylib:$$PYTHONPATH $(sphinx-build) -b man sphinx/ man/
	@echo "Documentation available in $$(pwd)/man/bout.1"

# Run doxygen, ignore if it fails (leading '-')
doxygen:
	-cd doxygen && doxygen Doxyfile

# Run breathe-apidoc, ignore if it fails (leading '-')
breathe-autogen: doxygen
	-breathe-apidoc -f -o sphinx/_breathe_autogen doxygen/bout/xml

clean:
	@echo "Cleaning up..."
	@$(RM) -f *.pdf
	@$(RM) -rf html/
	@$(RM) -rf build/
	@$(RM) -rf man/
	@$(RM) -rf sphinx/_breathe_autogen/
