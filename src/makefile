
BOUT_TOP	= ..

DIRS			= field fileio invert mesh physics precon solver sys facets

SOURCEC		= bout++.cxx
SOURCEH		= bout.hxx
CFLAGS		= -DMD5SUM=$(CHECKSUM) -DREVISION="$(REVISION)"
TARGET		= lib

# Checksum passed to bout++.cxx
CHECKSUM := $(shell cat $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.cxx $(DIR)/*.hxx)) | md5sum | head -c 32)

# Git revision passed to bout++.cxx
REVISION := $(shell git rev-parse HEAD 2> /dev/null || hg id -n 2> /dev/null || { echo "Unknown"; } )

# This is a little obscure, so here is an explanation:
# the first 'all: initial_message' sets the prerequisite 'inital_message'
# to be the first one called

all: initial_message

# Then the config file is included, which has all the other prerequisites for 'all'
include $(BOUT_TOP)/make.config

# Then set the last two preqrequisites for 'all'
all: lib end_message

checksum:
	@echo $(CHECKSUM)

# This to make sure build time is always printed correctly
initial_message:
	@touch $(SOURCEC)
	@echo ""
	@echo "----- Compiling BOUT++ -----" 
	@echo "CXX      = " $(CXX)
	@echo "CFLAGS   = " $(BOUT_FLAGS)
	@echo "CHECKSUM = " $(CHECKSUM)
	@echo "INCLUDE  = " $(BOUT_INCLUDE)
	@rm -f $(BOUT_TOP)/lib/libbout++.a

end_message:
	@echo "Finished"
	@echo ""
