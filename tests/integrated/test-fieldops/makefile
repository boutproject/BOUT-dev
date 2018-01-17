BOUT_TOP	= ../../..
SOURCEC = test.cxx
include $(BOUT_TOP)/make.config

test.o: fieldops.cxx

fieldops.cxx: gen_fieldops.py
	python3 $< > $@
	@#touch test.cxx

clean:: clean_log

clean_log:
	@rm run.log* -f
	@rm fieldops.cxx -f
	@rm test_fieldops.cxx -f
