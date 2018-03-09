BOUT_TOP=..

all:flexible.hxx

flexible.hxx:flexible.hxx.in.py
	PYTHONPATH+=$(BOUT_TOP)/src/field:$$PYTHONPATH ./$< > $@.tmp
	mv $@.tmp $@
