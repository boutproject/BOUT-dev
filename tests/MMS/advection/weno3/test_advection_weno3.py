#!/usr/bin/env python3

# requires: all_tests
# Requires: not make

import parent

parent.make_parent()

parent.nxlist = [64, 128]

parent.run_mms([["method=0 mesh:ddx:upwind=w3 mesh:ddz:upwind=w3", 3]])
