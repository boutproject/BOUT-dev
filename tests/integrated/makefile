# Makefile for example codes
# Works in the same way as the root Makefile

BOUT_TOP = ../..

LIST = test_suite_list

BUILD = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite_make --get-list)

CHECKING = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite --get-list)


include $(BOUT_TOP)/make.config

CHECKING_ = $(CHECKING:%=%_checking)

.PHONY: $(BUILD) $(CHECKING_) check

all: $(BUILD)

buildncheck: all check

$(BUILD):
	@$(MAKE) --no-print-directory -C $@

check: $(CHECKING_)
	@echo $(CHECKING_)

$(CHECKING_):
	@$(MAKE) --no-print-directory -C $(@:%_checking=%) runtest
