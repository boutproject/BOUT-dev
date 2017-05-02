# Makefile for example codes
# Works in the same way as the root Makefile

BOUT_TOP	= ..

list=test_suite_list

BUILD		= $(shell grep '^[^\#!]' $(list) )

CHECKING	= $(shell grep '^[^\#]'  $(list) |sed s/\!// )

include $(BOUT_TOP)/make.config


CHECKING_=$(CHECKING:%=%_checking)

.PHONY: $(BUILD) $(CHECKING_) check

all: $(BUILD)

buildncheck: all check

$(BUILD):
	@$(MAKE) --no-print-directory -C $@

check: $(CHECKING_)
	@echo $(CHECKING_)

$(CHECKING_):
	@$(MAKE) --no-print-directory -C $(@:%_checking=%) runtest
