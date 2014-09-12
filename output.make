
BOUT_TOP=$(PWD)

include make.config

.PHONY: cflags
cflags:
	@echo $(BOUT_INCLUDE) $(BOUT_FLAGS)

.PHONY: ldflags
ldflags:
	@echo $(LDFLAGS) $(BOUT_LIBS)

