BOUT_TOP = .

include make.config

srcdirs = src


all: make.config $(srcdirs)

.PHONY:$(srcdirs)
$(srcdirs):
	@$(MAKE) --no-print-directory -C $@

config:
	./configure

clean::
	@for pp in $(srcdirs); do echo Cleaning $$pp; $(MAKE) --no-print-directory -C $$pp clean; cd; done
	@echo Cleaning examples
	@$(MAKE) --no-print-directory -C examples clean;

.PHONY: examples
examples: all
	@$(MAKE) -C examples;
