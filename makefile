BOUT_TOP  = $(PWD)

DIRS      = src

DIRS_CLEAN =  tests/integrated tests/unit tests/MMS

TARGET   ?= libfast

include make.config

# For compatibility ignore shared
shared: libfast

######################################################################
# Tests
######################################################################

check-unit-tests: libfast
	@export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH}; $(MAKE) --no-print-directory -C tests/unit check

check-mms-tests: libfast
	+@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} \
		OMPI_MCA_rmaps_base_oversubscribe=yes ./test_suite_make
	+@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} \
		OMPI_MCA_rmaps_base_oversubscribe=yes ./test_suite

check-mms-tests-all: libfast
	+@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} \
		OMPI_MCA_rmaps_base_oversubscribe=yes ./test_suite_make --all
	+@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} \
		OMPI_MCA_rmaps_base_oversubscribe=yes ./test_suite --all

check-integrated-tests: libfast
	+@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} \
		OMPI_MCA_rmaps_base_oversubscribe=yes ./test_suite_make
	+@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} \
		OMPI_MCA_rmaps_base_oversubscribe=yes ./test_suite

check-integrated-tests-all: libfast
	+@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} \
		OMPI_MCA_rmaps_base_oversubscribe=yes ./test_suite_make --all
	+@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} \
		OMPI_MCA_rmaps_base_oversubscribe=yes ./test_suite --all


check: check-unit-tests check-integrated-tests check-mms-tests
check-all: check-unit-tests check-integrated-tests-all check-mms-tests-all

build-check-unit-tests: libfast
	@$(MAKE) --no-print-directory -C tests/unit

build-check-mms-tests: libfast
	$(MAKE) --no-print-directory -C tests/MMS
build-check-mms-tests-all: libfast
	$(MAKE) --no-print-directory -C tests/MMS allall

build-check-integrated-tests: libfast
	$(MAKE) --no-print-directory -C tests/integrated
build-check-integrated-tests-all: libfast
	$(MAKE) --no-print-directory -C tests/integrated allall


build-check: build-check-integrated-tests build-check-mms-tests build-check-unit-tests
build-check-all: build-check-integrated-tests-all build-check-mms-tests-all build-check-unit-tests

# Build the .mo files needed for Natural Language Support (gettext)
.PHONY: locale
locale:
	$(MAKE) -C locale

######################################################################
# Releases
######################################################################

.PHONY: dist changelog

# Makes the tarball BOUT++-v<version>.tar.gz
dist:
	@bin/bout-archive-helper.sh v$(BOUT_VERSION)

CHANGELOG_ERR_MESSAGE := "Run like: make changelog TOKEN=<token> LAST_VERSION=vX.Y.Z RELEASE_BRANCH=master|next"

# Updates CHANGELOG.md, needs some arguments:
#
#     make changelog TOKEN=<token> LAST_VERSION=vX.Y.Z RELEASE_BRANCH=master|next
#
# Note: You should probably only run this if you are a maintainer (and
# also know what you're doing)!
changelog:
ifndef TOKEN
	$(error $(CHANGELOG_ERR_MESSAGE))
endif
ifndef LAST_VERSION
	$(error $(CHANGELOG_ERR_MESSAGE))
endif
ifndef RELEASE_BRANCH
	$(error $(CHANGELOG_ERR_MESSAGE))
endif
	github_changelog_generator -t $(TOKEN) --since-tag \
        $(LAST_VERSION) --no-issues --max-issues 700 \
        --base CHANGELOG.md --future-release v$(BOUT_VERSION) \
        --release-branch $(RELEASE_BRANCH) \
        --user boutproject --project BOUT-dev
