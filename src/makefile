
BOUT_TOP	= ..

DIRS			= field fileio invert mesh physics precon solver sys

SOURCEC		= bout++.cxx
SOURCEH		= bout.hxx
CXXFLAGS       += -DMD5SUM=$(CHECKSUM) -DREVISION="$(REVISION)"
ifndef TARGET
TARGET		= lib
endif

# Checksum passed to bout++.cxx
CHECKSUM := $(shell cat $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.cxx $(DIR)/*.hxx)) | md5sum | head -c 32 )

# Git revision passed to bout++.cxx
REVISION := $(shell git rev-parse HEAD 2> /dev/null || hg id -n 2> /dev/null || { echo "Unknown"; } )

# This is a little obscure, so here is an explanation:
# the first 'all: initial_message' sets the prerequisite 'inital_message'
# to be the first one called

all: initial_message

# Then the config file is included, which has all the other prerequisites for 'all'
include $(BOUT_TOP)/make.config


.PHONY: bout++.cxx

libfast:
	@echo "Recreating libbout++.a"
	@rm -f $(BOUT_TOP)/lib/libbout++.a
	@$(AR) $(ARFLAGS) $(LIB) $(shell find $(BOUT_TOP)/src -name \*.o -type f -print 2> /dev/null)
	@#$(RANLIB) $(LIB)

# Then set the last two preqrequisites for 'all'
all: lib end_message

checksum:
	@echo $(CHECKSUM)

# This to make sure build time is always printed correctly
.PHONY: bout++.cxx

initial_message:
	@echo ""
	@echo "----- Compiling BOUT++ -----" 
	@echo "CXX      = " $(CXX)
	@echo "FLAGS    = " $(BOUT_FLAGS)
	@echo "CHECKSUM = " $(CHECKSUM)
	@echo "INCLUDE  = " $(BOUT_INCLUDE)
	@rm -f $(BOUT_TOP)/lib/libbout++.a

end_message:
	@echo "Finished"
	@echo ""
