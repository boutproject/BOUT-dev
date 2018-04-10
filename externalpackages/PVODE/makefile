# Makefile for example codes
# Works in the same way as the root Makefile

TARGET?=all

DIRS		= source precon

all: $(DIRS)

.PHONY: $(DIRS)

$(DIRS):
	@$(MAKE) -s --no-print-directory TARGET=$(TARGET) -C $@ $(TARGET)

clean::
	@for pp in $(DIRS); do echo "  " $$pp cleaned; $(MAKE) --no-print-directory -C $$pp clean; done
