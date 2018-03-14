BOUT_TOP=..

all:flexible.hxx

flexible.hxx:flexible.hxx.in.py
	@echo "  Generating $@"
	@PYTHONPATH+=$(BOUT_TOP)/src/field:$$PYTHONPATH ./$< > $@.tmp || \
               (ex=$? ; echo "  touch $@ to ignore error and use old version" ; exit $ex)
	@mv $@.tmp $@
