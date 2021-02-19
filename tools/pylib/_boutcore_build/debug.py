#!/usr/bin/env python3

print("Trying import")
import debug
import sys

print(sys.argv)
print("Trying init")
debug.blas()
print("Trying to get rank")
print(debug.rank())
print("Trying finalize")
debug.finish()
