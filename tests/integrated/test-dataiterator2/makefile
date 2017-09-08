
BOUT_TOP	= ../../..

SOURCEC		= test_dataiterator.cxx

include $(BOUT_TOP)/make.config

test3d.cxx:gen.py
	python2 $< > $@

test_dataiterator.o:test3d.cxx

BOUT_FLAGS += -g -O0  
