
BOUT_TOP	= ..

DIRS			= field fileio invert mesh physics solver sys cbout

SOURCEC		= bout++.cxx
SOURCEH		= bout.hxx
CXXFLAGS       += -DMD5SUM=$(CHECKSUM) -DREVISION="$(REVISION)"
TARGET	       ?= lib

# Checksum passed to bout++.cxx
CHECKSUM := $(shell cat $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.cxx $(DIR)/*.hxx)) | md5sum | head -c 32 )

# Git revision passed to bout++.cxx
REVISION := $(shell git rev-parse HEAD 2> /dev/null || hg id -n 2> /dev/null || { echo "Unknown"; } )

# From the legacy build script:
# This is a little obscure, so here is an explanation:
# the first 'all: initial_message' sets the prerequisite 'inital_message'
# to be the first one called
all: initial_message

# Then the config file is included, which has all the other
# prerequisites for 'all'
include $(BOUT_TOP)/make.config

# dummy target - must not be PHONY, to prevent unneeded rebuild of bout++.o
# A target is needed to make sure that the directories are build first
# It must be a real and existing target, otherwise the depending
# targets are always updated as the requirement (.dummy) does not
# exist.
.dummy:$(DIRS)
	@#empty recipe is needed

# Rebuild bout++ to get the correct time of the build - but only if
# any other file was rebuild
# .dummy dependency to make sure the library is not build before the
# directories are finished and further to make sure we actually know
# whether the file needs to be rebuild - which we can only know for
# certain if all directories are finished
bout++.o: $(BOUT_TOP)/lib/.last.o.file .dummy

# The recipie could be removed, as it only catches the case of
# out-of-date make.config
# The rest is needed to tell make that .last.o.file will be created
$(BOUT_TOP)/lib/.last.o.file: .dummy
	@test -f $@ || \
	   echo "make.config is out of date - run ./configure"
	@test -f $@ || \
	   exit 1

# make .libfast a real target - similar to .dummy to avoid that a
# rebuild is forced.
libfast: .libfast

.libfast: bout++.o
	@echo "Recreating libbout++.a"
	@rm -f $(BOUT_TOP)/lib/libbout++.a
	@$(AR) $(ARFLAGS) $(LIB) $(shell find $(BOUT_TOP)/src -name \*.o -type f -print 2> /dev/null)
	@touch .libfast
	@#$(RANLIB) $(LIB)

# From the legacy build script:
# Then set the last two preqrequisites for 'all'
all: lib end_message

# Only needed for legacy build script
ifneq ("$(TARGET)","libfast")
# This to make sure build time is always printed correctly
.PHONY: bout++.cxx
endif

checksum:
	@echo $(CHECKSUM)

# From the legacy build script
initial_message:
	@echo ""
	@echo "----- Compiling BOUT++ -----" 
	@echo "CXX      = " $(CXX)
	@echo "FLAGS    = " $(BOUT_FLAGS)
	@echo "CHECKSUM = " $(CHECKSUM)
	@echo "INCLUDE  = " $(BOUT_INCLUDE)
	@rm -f $(BOUT_TOP)/lib/libbout++.a

# From the legacy build script
end_message:
	@echo "Finished"
	@echo ""
