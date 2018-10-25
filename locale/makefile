######################################################################
# Internationalisation
######################################################################

# List of supported locales
LANGUAGES = fr zh_TW zh_CN

# This target is the most commonly used, to build all the .mo files
.PHONY: locale
locale: $(LANGUAGES:%=%/LC_MESSAGES/libbout.mo)

# Update all .po files. Note that they will then need manually
# editing to fill in (or delete) missing translations.
.PHONY: po-all
all-po: $(LANGUAGES:%=%/libbout.po)

# Create the template file, combining all source files
libbout.pot: FORCE
	xgettext --keyword=_ --language=c++ --add-comments --sort-output -o $@ `find ../ -name "*.[ch]xx"`

# Update a .po file
%/libbout.po: libbout.pot
	msgmerge --update $@ $<

# Update a .mo file, ensuring that the directory exists
# Note: Because the .po files must be updated manually,
#       don't automatically generate the .po files
%/LC_MESSAGES/libbout.mo: FORCE
	@mkdir -p $(@D)
	msgfmt --output-file=$@ $*/libbout.po

# Shortcut target for building single language
# e.g. "make po-fr" 
# Note: invoking make since otherwise the intermediate file is deleted
po-%:
	$(MAKE) $*/libbout.po
mo-%:
	$(MAKE) $*/LC_MESSAGES/libbout.mo

# Dummy target to force rules to be run
FORCE:
