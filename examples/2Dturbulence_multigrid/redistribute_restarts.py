#!/usr/bin/env python

from boutdata import restart
from sys import argv, exit

npes = int(argv[1])

try:
    restarts_directory = argv[2]
except IndexError:
    restarts_directory = "restarts_256x256"

restart.redistribute(npes, path=restarts_directory, output="data", myg=0)

exit(0)
