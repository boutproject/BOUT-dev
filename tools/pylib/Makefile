all: python2 python3
python%:
	$(MAKE) -C _boutcore_build clean
	PY=$@ $(MAKE) -C _boutcore_build pseudoinstall
