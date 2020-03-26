######################################################################
# Internationalisation
######################################################################

# Find all .po files and work out the corresponding .mo file name
PO_FILES = $(shell ls */libbout.po)
MO_FILES = $(PO_FILES:%/libbout.po=%/LC_MESSAGES/libbout.mo)

# This target is the most commonly used, to build all the .mo files
.PHONY: locale
locale: $(MO_FILES)

# Create the template file, combining all source files
libbout.pot: FORCE
	xgettext --keyword=_ --language=c++ --add-comments --sort-output -o $@ `find ../ -name "[^.#]*.[ch]xx"`

# Update a .po file
# If it doesn't exist then create; if it exists then update
# Note that it will need manually editing to translate or delete new entries
%/libbout.po: libbout.pot
	@mkdir -p $(@D)
	@if [ -a $@ ]; then echo "Updating " $@; msgmerge --update $@ $< ; else echo "Creating " $@; msginit --input=libbout.pot --locale=$* --output=$@; fi;
	@echo "Now edit " $@ " to translate or delete new entries"
	@echo "then run 'make' to generate the .mo file"

# Update a .mo file, ensuring that the directory exists
# Note: Because the .po files must be updated manually,
#       don't automatically generate the .po files
%/LC_MESSAGES/libbout.mo: FORCE
	@mkdir -p $(@D)
	@echo "Building language: " $*
	@msgfmt --output-file=$@ $*/libbout.po

# Shortcut target for building single language
# e.g. "make locale-fr" 
# Note: invoking make since otherwise the intermediate file is deleted
locale-%:
	$(MAKE) $*/libbout.po
mo-%:
	$(MAKE) $*/LC_MESSAGES/libbout.mo

# Dummy target to force rules to be run
FORCE:

update-all: libbout.pot
	for d in */ ; do $(MAKE) --no-print-directory $${d}libbout.po; done
	for d in */ ; do $(MAKE) --no-print-directory $${d}LC_MESSAGES/libbout.mo; done
	grep Last-Translator */*.po
