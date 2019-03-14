
include make.config

# These expand to literals of the same name, which are then
# set in the bout.config script. This enables
# the final paths to the include and lib files to be set
# in the make/make install targets

BOUT_INCLUDE_PATH="\$$BOUT_INCLUDE_PATH"
BOUT_LIB_PATH="\$$BOUT_LIB_PATH"
MPARK_VARIANT_INCLUDE_PATH="\$$MPARK_VARIANT_INCLUDE_PATH"

.PHONY: cflags
cflags:
	@echo $(BOUT_INCLUDE) $(BOUT_FLAGS)

.PHONY: ldflags
ldflags:
	@echo $(LDFLAGS) $(BOUT_LIBS)

