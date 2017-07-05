all: python2 python3
python%:
	make -C _boutcore_build clean
	PY=$@ make -C _boutcore_build pseudoinstall
