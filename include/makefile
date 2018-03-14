BOUT_TOP=..

all: derivs.hxx flexible.hxx

derivs.hxx: derivs.hxx.in.py derivs.hxx.in.jinja
# run first dependency and write output to target
# on failure delete output and report error
	@echo "  Generating $@"
	@./$< > $@.tmp || (fail=$$?; echo "touch $@ to ignore failed generation" ; exit $$fail)
	@mv $@.tmp $@
	@clang-format -i $@ || echo "Formatting failed"

flexible.hxx:flexible.hxx.in.py
	@echo "  Generating $@"
	@PYTHONPATH+=$(BOUT_TOP)/src/field:$$PYTHONPATH ./$< > $@.tmp || \
               (ex=$? ; echo "  touch $@ to ignore error and use old version" ; exit $ex)
	@clang-format -i $@.tmp || echo "Warning: Formatting failed - please do not add file to repository!"
	@mv $@.tmp $@
